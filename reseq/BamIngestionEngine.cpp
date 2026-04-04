#include "BamIngestionEngine.h"
using reseq::BamIngestionEngine;

#include "DataStats.h"
using reseq::DataStats;
using reseq::Vect;

#include <algorithm>
using std::max;
using std::min;
#include <chrono>
// include <array>
using std::array;
#include <cstring>
using std::strlen;
#include <mutex>
using std::lock_guard;
using std::mutex;
using std::unique_lock;
// include <string>
using std::string;
// include <utility>
using std::pair;
#include <thread>
// include <vector>
using std::vector;

#include "CMakeConfig.h"
#include "logging.hpp"

// include <seqan/bam_io.h>
using seqan::atEnd;
using seqan::BamAlignmentRecord;
using seqan::BamFileIn;
using seqan::BamHeader;
using seqan::contigNames;
using seqan::Dna;
using seqan::Dna5;
using seqan::Exception;
using seqan::readHeader;
using seqan::readRecord;

// include "utilities.hpp"
using reseq::utilities::at;
using reseq::utilities::Complement;
using reseq::utilities::ConstIupacStringReverseComplement;
using reseq::utilities::Divide;
using reseq::utilities::MeanWithRoundingToFirst;
using reseq::utilities::Percent;
using reseq::utilities::ReversedConstCharString;
using reseq::utilities::ReversedConstCigarString;
using reseq::utilities::SafePercent;
using reseq::utilities::SetToMax;
using reseq::utilities::SetToMin;

uint64_t BamIngestionEngine::RecordHasher::operator()(CoverageStats::FullRecord* const& record) const {
    // FNV-1a hash: https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
    const auto cdata = reinterpret_cast<unsigned char*>(toCString(record->record_.qName));
    uint64_t hash = 14695981039346656037ull;
    for (auto i = length(record->record_.qName); i--;) {
        const auto next = std::size_t{cdata[i]};
        hash = (hash ^ next) * 1099511628211ull;
    }

    return hash;
}

bool BamIngestionEngine::RecordEqual::operator()(CoverageStats::FullRecord* const& lhs,
                                                 CoverageStats::FullRecord* const& rhs) const {
    return lhs->record_.qName == rhs->record_.qName;
}

inline bool BamIngestionEngine::PotentiallyValidGeneral(const seqan::BamAlignmentRecord& record,
                                                        const DataStats& target) const {
    return !hasFlagNextUnmapped(record) && !hasFlagUnmapped(record) && record.rNextId == record.rID &&
           !target.reference_->ReferenceSequenceExcluded(record.rID);
}
inline bool BamIngestionEngine::PotentiallyValidFirst(const seqan::BamAlignmentRecord& record_first,
                                                      const DataStats& target) const {
    return PotentiallyValidGeneral(record_first, target) &&
           record_first.beginPos + target.maximum_insert_length_ >= record_first.pNext;
}
inline bool BamIngestionEngine::PotentiallyValidSecond(const seqan::BamAlignmentRecord& record_second,
                                                       const DataStats& target) const {
    return PotentiallyValidGeneral(record_second, target) &&
           record_second.pNext + target.maximum_insert_length_ >= record_second.beginPos;
}
inline bool BamIngestionEngine::PotentiallyValid(const seqan::BamAlignmentRecord& record,
                                                 const DataStats& target) const {
    if (record.beginPos < record.pNext) {
        return PotentiallyValidFirst(record, target);
    } else {
        return PotentiallyValidSecond(record, target);
    }
}

bool BamIngestionEngine::IsSecondRead(CoverageStats::FullRecord* record, CoverageStats::FullRecord*& record_first,
                                      CoverageStats::CoverageBlock*& block, DataStats& target) {
    auto insert_info = first_read_records_.insert(record);
    if (insert_info.second) {
        // Insertion worked, so it is the first read
        if (PotentiallyValidFirst(record->record_, target)) {
            // Use record->record_.beginPos-maximum_read_length_on_reference_ as start, because of potentially
            // soft-clipped bases at the beginning
            target.coverage_.AddFragment(record->record_.rID,
                                         (record->record_.beginPos > target.maximum_read_length_on_reference_
                                              ? record->record_.beginPos - target.maximum_read_length_on_reference_
                                              : 0),
                                         block);
        }

        return false;
    } else {
        // Insertion failed, so it is the second read
        record_first = *insert_info.first;
        first_read_records_.erase(insert_info.first);

        return true;
    }
}

bool BamIngestionEngine::CheckForAdapters(const seqan::BamAlignmentRecord& record_first,
                                          const seqan::BamAlignmentRecord& record_second, DataStats& target) {
    uintReadLen adapter_position_first, adapter_position_second;
    bool adapter_detected = target.adapters_.Detect(adapter_position_first, adapter_position_second, record_first,
                                                    record_second, *target.reference_);

    if (adapter_detected) {
        auto insert_length = MeanWithRoundingToFirst(adapter_position_first, adapter_position_second);
        target.fragment_distribution_.AddInsertLengths(insert_length);
        ++target.tmp_non_mapped_read_lengths_by_fragment_length_.at(0).at(insert_length).at(length(record_first.seq));
        ++target.tmp_non_mapped_read_lengths_by_fragment_length_.at(1).at(insert_length).at(length(record_second.seq));

        return true;
    } else {
        return false;
    }
}

bool BamIngestionEngine::EvalReferenceStatistics(CoverageStats::FullRecord* record, uintTempSeq template_segment,
                                                 CoverageStats::CoverageBlock* coverage_block, DataStats& target) {
    // Get gc on reference
    auto cov_block = coverage_block; // coverage_block might be changed in the next step so keep the original for later
    uintSeqLen coverage_pos = target.coverage_.GetStartPos(record->from_ref_pos_, cov_block);
    array<uintNucCount, 5> seq_content_reference = {0, 0, 0, 0, 0};
    for (auto ref_pos = record->from_ref_pos_; ref_pos < record->to_ref_pos_; ++ref_pos) {
        if (cov_block->coverage_.at(ref_pos - cov_block->start_pos_).valid_) {
            ++seq_content_reference.at(at(target.reference_->ReferenceSequence(record->record_.rID), ref_pos));
        } else {
            ++seq_content_reference.at(4);
        }
        target.coverage_.IncrementPos(coverage_pos, cov_block);
    }
    auto gc_percent = SafePercent(seq_content_reference.at(1) + seq_content_reference.at(2),
                                  record->to_ref_pos_ - record->from_ref_pos_ - seq_content_reference.at(4));
    record->reference_gc_ = gc_percent;

    uintSeqLen ref_pos;
    if (hasFlagRC(record->record_)) {
        coverage_pos = target.coverage_.GetStartPos(record->to_ref_pos_ - 1, coverage_block);
        ref_pos = record->to_ref_pos_ - 1;
    } else {
        coverage_pos = target.coverage_.GetStartPos(record->from_ref_pos_, coverage_block);
        ref_pos = record->from_ref_pos_;
    }

    ConstIupacStringReverseComplement reversed_seq(record->record_.seq);
    ReversedConstCharString reversed_qual(record->record_.qual);
    ReversedConstCigarString rev_cigar(record->record_.cigar);

    uintReadLen read_pos(0), read_pos_ref(0);
    Dna ref_base;
    Dna5 base;
    uintBaseCall last_base(5);
    uintQual qual;

    array<uintNucCount, 5> seq_content_mapped = {0, 0, 0, 0, 0};
    uintReadLen indel_pos(0);
    uintInDelType indel_type(0);

    for (const auto& cigar_element : (hasFlagRC(record->record_) ? rev_cigar : record->record_.cigar)) {
        switch (cigar_element.operation) {
        case 'M':
        case '=':
        case 'X':
        case 'S': // Treat soft-clipping as match, so that bwa and bowtie2 behave the same
            for (auto i = cigar_element.count; i--;) {
                if (record->from_ref_pos_ <= ref_pos && ref_pos < record->to_ref_pos_) {
                    // We are currently not comparing the adapter to the reference
                    if (hasFlagRC(record->record_)) {
                        base = at(reversed_seq, read_pos);
                        qual = at(reversed_qual, read_pos) - target.phred_quality_offset_;
                        ref_base =
                            Complement::Dna5(at(target.reference_->ReferenceSequence(record->record_.rID), ref_pos));
                    } else {
                        base = at(record->record_.seq, read_pos);
                        qual = at(record->record_.qual, read_pos) - target.phred_quality_offset_;
                        ref_base = at(target.reference_->ReferenceSequence(record->record_.rID), ref_pos);
                    }

                    if (!coverage_block->coverage_.at(coverage_pos).valid_) {
                        ++seq_content_mapped.at(4);
                    } else {
                        ++target.tmp_sequence_content_reference_.at(template_segment)
                              .at(hasFlagRC(record->record_))
                              .at(ref_base)
                              .at(read_pos_ref);

                        target.errors_.AddBasePlotting(template_segment, ref_base, base, qual, last_base);
                        target.errors_.AddInDel(indel_type, last_base, ErrorStats::InDelDef::kNoInDel, indel_pos,
                                                read_pos, gc_percent);

                        ++seq_content_mapped.at(base);
                    }

                    last_base = base;

                    ++read_pos;
                    ++read_pos_ref;

                    if (hasFlagRC(record->record_)) {
                        target.coverage_.AddReverse(coverage_pos, coverage_block, base);

                        target.coverage_.DecrementPos(coverage_pos, coverage_block);
                        --ref_pos;
                    } else {
                        target.coverage_.AddForward(coverage_pos, coverage_block, base);

                        target.coverage_.IncrementPos(coverage_pos, coverage_block);
                        ++ref_pos;
                    }

                    indel_type = 0;
                    indel_pos = 0;
                }
            }
            break;
        case 'N':
        case 'D':
            for (indel_pos = 0; indel_pos < cigar_element.count; ++indel_pos) {
                if (record->from_ref_pos_ <= ref_pos && ref_pos < record->to_ref_pos_) {
                    // We are currently not comparing the adapter to the reference
                    if (hasFlagRC(record->record_)) {
                        ref_base =
                            Complement::Dna5(at(target.reference_->ReferenceSequence(record->record_.rID), ref_pos));
                    } else {
                        ref_base = at(target.reference_->ReferenceSequence(record->record_.rID), ref_pos);
                    }

                    if (coverage_block->coverage_.at(coverage_pos).valid_) {
                        ++target.tmp_sequence_content_reference_.at(template_segment)
                              .at(hasFlagRC(record->record_))
                              .at(ref_base)
                              .at(read_pos_ref);

                        target.errors_.AddInDel(indel_type, last_base, ErrorStats::InDelDef::kDeletion, indel_pos,
                                                read_pos, gc_percent);
                    }

                    ++read_pos_ref;

                    if (hasFlagRC(record->record_)) {
                        target.coverage_.DecrementPos(coverage_pos, coverage_block);
                        --ref_pos;
                    } else {
                        target.coverage_.IncrementPos(coverage_pos, coverage_block);
                        ++ref_pos;
                    }
                    indel_type = 1;
                }
            }
            break;
        case 'I':
            if (record->from_ref_pos_ <= ref_pos && ref_pos < record->to_ref_pos_) {
                if (hasFlagRC(record->record_)) {
                    ref_base = Complement::Dna5(at(target.reference_->ReferenceSequence(record->record_.rID), ref_pos));
                } else {
                    ref_base = at(target.reference_->ReferenceSequence(record->record_.rID), ref_pos);
                }

                if (!coverage_block->coverage_.at(coverage_pos).valid_) {
                    read_pos += cigar_element.count;
                    indel_type = 0;
                } else {
                    for (indel_pos = 0; indel_pos < cigar_element.count; ++indel_pos) {
                        if (hasFlagRC(record->record_)) {
                            base = at(reversed_seq, read_pos);
                        } else {
                            base = at(record->record_.seq, read_pos);
                        }

                        target.errors_.AddInDel(indel_type, last_base,
                                                static_cast<ErrorStats::InDelDef>(static_cast<uintBaseCall>(base) + 2),
                                                indel_pos, read_pos, gc_percent);
                        ++read_pos;
                        indel_type = 0;
                    }
                }
                seq_content_mapped.at(4) += cigar_element.count; // Ignore inserted bases for gc_percent
            }
            break;
        default:
            break;
        }
    }

    ++target.tmp_gc_read_content_reference_.at(template_segment).at(gc_percent);
    ++target.tmp_gc_read_content_mapped_.at(template_segment)
          .at(SafePercent(seq_content_mapped.at(1) + seq_content_mapped.at(2),
                          length(record->record_.seq) - seq_content_mapped.at(4)));

    return true;
}

bool BamIngestionEngine::EvalRecord(pair<CoverageStats::FullRecord*, CoverageStats::FullRecord*> record,
                                    ThreadData& thread, DataStats& target) {
    // Start handling pair by determining tile(and counting it) and template segment
    uintTileId tile_id;
    if (!target.tiles_.GetTileId(tile_id, record.first->record_.qName)) {
        return false;
    }

    uintTempSeq template_segment(0);
    if (hasFlagLast(record.second->record_)) {
        template_segment = 1;
    }

    uintTempSeq strand;
    if (hasFlagUnmapped(record.first->record_)) {
        if (hasFlagUnmapped(record.second->record_)) {
            strand = 2; // Invalid strand as we cannot know it without mapping
        } else {
            strand = (hasFlagRC(record.second->record_) && hasFlagFirst(record.second->record_)) ||
                     (!hasFlagRC(record.second->record_) && hasFlagLast(record.second->record_));
        }
    } else {
        strand = (hasFlagRC(record.first->record_) && hasFlagFirst(record.first->record_)) ||
                 (!hasFlagRC(record.first->record_) && hasFlagLast(record.first->record_));
    }

    // Base level stats
    uintQual paired_seq_qual(0);
    target.read_sequence_stats_.EvalBaseLevelStats(record.first, (template_segment + 1) % 2, strand, tile_id,
                                                   paired_seq_qual, target.qualities_, target.phred_quality_offset_,
                                                   target.maximum_quality_);
    target.read_sequence_stats_.EvalBaseLevelStats(record.second, template_segment, strand, tile_id, paired_seq_qual,
                                                   target.qualities_, target.phred_quality_offset_,
                                                   target.maximum_quality_);

    // Mapping and Reference stats
    if (hasFlagUnmapped(record.first->record_) || hasFlagNextUnmapped(record.first->record_)) {
        if (!hasFlagUnmapped(record.first->record_)) {
            target.read_sequence_stats_.IncrementMappingQuality(2, record.first->record_.mapQ);
        } else if (!hasFlagUnmapped(record.second->record_)) {
            target.read_sequence_stats_.IncrementMappingQuality(2, record.second->record_.mapQ);
        }

        if (CheckForAdapters(record.first->record_, record.second->record_, target)) {
            reads_in_unmapped_pairs_with_adapters_ += 2;
        } else {
            reads_in_unmapped_pairs_without_adapters_ += 2;
        }
    } else {
        if (!hasFlagAllProper(record.first->record_)) {
            target.read_sequence_stats_.IncrementMappingQuality(1, record.first->record_.mapQ);
            target.read_sequence_stats_.IncrementMappingQuality(1, record.second->record_.mapQ);
        } else {
            target.read_sequence_stats_.IncrementMappingQuality(0, record.first->record_.mapQ);
            target.read_sequence_stats_.IncrementMappingQuality(0, record.second->record_.mapQ);
        }

        if (!target.QualitySufficient(record.first->record_) || !target.QualitySufficient(record.second->record_)) {
            if (target.reference_->ReferenceSequenceExcluded(record.first->record_.rID) ||
                target.reference_->ReferenceSequenceExcluded(record.second->record_.rID)) {
                reads_on_too_short_fragments_ += 2;
            } else {
                uintSeqLen start_pos_first, end_pos_first, start_pos_second, end_pos_second;
                target.GetReadPosOnReference(start_pos_first, end_pos_first, record.first->record_);
                target.GetReadPosOnReference(start_pos_second, end_pos_second, record.second->record_);
                if (target.reference_->FragmentExcluded(thread.last_exclusion_region_id_,
                                                        thread.last_exclusion_ref_seq_, record.first->record_.rID,
                                                        start_pos_first, end_pos_first) ||
                    target.reference_->FragmentExcluded(thread.last_exclusion_region_id_,
                                                        thread.last_exclusion_ref_seq_, record.second->record_.rID,
                                                        start_pos_second, end_pos_second)) {
                    reads_in_excluded_regions_ += 2;
                } else {
                    if (CheckForAdapters(record.first->record_, record.second->record_,
                                         target)) { // We don't trust the mapping, so look for adapters
                                                    // like in unmapped case
                        reads_with_low_quality_with_adapters_ += 2;
                    } else {
                        reads_with_low_quality_without_adapters_ += 2;
                    }
                }
            }
        } else {
            if (record.first->record_.rID == record.second->record_.rID) {
                if (target.reference_->ReferenceSequenceExcluded(record.first->record_.rID)) {
                    reads_on_too_short_fragments_ += 2;
                } else {
                    uintSeqLen start_pos_first, start_pos_second, end_pos_first, end_pos_second;
                    target.GetReadPosOnReference(start_pos_first, end_pos_first, record.first->record_);
                    target.GetReadPosOnReference(start_pos_second, end_pos_second, record.second->record_);

                    if (start_pos_second < start_pos_first) {
                        // Switch records as the soft-trimming changed order in bam sort
                        auto tmp = record.first;
                        record.first = record.second;
                        record.second = tmp;

                        template_segment = (template_segment + 1) % 2;

                        auto tmp2 = start_pos_first;
                        start_pos_first = start_pos_second;
                        start_pos_second = tmp2;
                        tmp2 = end_pos_first;
                        end_pos_first = end_pos_second;
                        end_pos_second = tmp2;
                    }

                    bool adapter_detected(false);
                    if (hasFlagRC(record.first->record_) && !hasFlagRC(record.second->record_) &&
                        end_pos_first > start_pos_second) {
                        adapter_detected = true;
                    }

                    if (target.InProperDirection(record.first->record_, end_pos_second, start_pos_first,
                                                 start_pos_second) ||
                        adapter_detected) {
                        if ((!adapter_detected && target.reference_->FragmentExcluded(
                                                      thread.last_exclusion_region_id_, thread.last_exclusion_ref_seq_,
                                                      record.first->record_.rID, start_pos_first, end_pos_second)) ||
                            (adapter_detected && target.reference_->FragmentExcluded(
                                                     thread.last_exclusion_region_id_, thread.last_exclusion_ref_seq_,
                                                     record.first->record_.rID, start_pos_second, end_pos_first))) {
                            reads_in_excluded_regions_ += 2;
                        } else {
                            reads_used_ += 2;

                            if (adapter_detected) {
                                uintReadLen adapter_position_first, adapter_position_second;
                                // Reads in proper direction for adapters (reverse-forward) and overlap is large enough
                                // (minimum half of the read length)
                                target.adapters_.Detect(adapter_position_first, adapter_position_second,
                                                        record.first->record_, record.second->record_,
                                                        *target.reference_,
                                                        true); // Call it just to insert the adapter stats, position is
                                                               // better known from mapping
                            }

                            CoverageStats::CoverageBlock* coverage_block;
                            // Pair stats
                            uintSeqLen insert_length;
                            if (adapter_detected) {
                                insert_length = end_pos_first - start_pos_second;

                                target.fragment_distribution_.AddGCContent(target.reference_->GCContent(
                                    record.second->record_.rID, start_pos_second, end_pos_first));
                                target.fragment_distribution_.FillInOutskirtContent(
                                    *target.reference_, record.second->record_, start_pos_second,
                                    end_pos_first); // Use second record here as it is forward, so if it is sequenced
                                                    // second, fragment is reversed
                                coverage_block =
                                    target.coverage_.FindBlock(record.second->record_.rID, start_pos_second);

                                record.first->from_ref_pos_ = start_pos_second;
                                record.second->from_ref_pos_ = start_pos_second;
                                record.first->to_ref_pos_ = end_pos_first;
                                record.second->to_ref_pos_ = end_pos_first;

                                target.fragment_distribution_.AddFragmentSite(
                                    record.first->record_.rID, insert_length, record.first->from_ref_pos_,
                                    template_segment, *target.reference_,
                                    thread.fragment_distribution_); // If second record(must be forward to be accepted
                                                                    // due to adapter) is sequenced second, fragment is
                                                                    // reversed
                            } else {
                                insert_length = end_pos_second - start_pos_first;

                                target.fragment_distribution_.AddGCContent(target.reference_->GCContent(
                                    record.first->record_.rID, start_pos_first, end_pos_second));
                                target.fragment_distribution_.FillInOutskirtContent(
                                    *target.reference_, record.first->record_, start_pos_first, end_pos_second);
                                coverage_block = target.coverage_.FindBlock(record.first->record_.rID, start_pos_first);

                                record.first->from_ref_pos_ = start_pos_first;
                                record.second->from_ref_pos_ = start_pos_second;
                                record.first->to_ref_pos_ = end_pos_first;
                                record.second->to_ref_pos_ = end_pos_second;

                                target.fragment_distribution_.AddFragmentSite(
                                    record.first->record_.rID, insert_length, record.first->from_ref_pos_,
                                    (template_segment + 1) % 2, *target.reference_,
                                    thread.fragment_distribution_); // If first record(must be forward to be accepted
                                                                    // without adapters) is sequenced second, fragment
                                                                    // is reversed
                            }

                            target.fragment_distribution_.AddAbundance(record.first->record_.rID);
                            target.fragment_distribution_.AddInsertLengths(
                                insert_length); // We do not just use the template length from field 9 of the bam file,
                                                // because bwa stores there the right-most position of the reverse read
                                                // - the left-most position of the forward read, which is not the
                                                // fragment length in the case of the reverse read being positioned
                                                // before the forward read. Bowtie2 gives the correct length here, but I
                                                // don't want to sacrifice compatibility so easily

                            ++target.tmp_read_lengths_by_fragment_length_.at(0)
                                  .at(insert_length)
                                  .at(length(record.first->record_.seq));
                            ++target.tmp_read_lengths_by_fragment_length_.at(1)
                                  .at(insert_length)
                                  .at(length(record.second->record_.seq));

                            record.first->tile_id_ = tile_id;
                            record.second->tile_id_ = tile_id;
                            record.first->fragment_length_ = insert_length;
                            record.second->fragment_length_ = insert_length;

                            if (!EvalReferenceStatistics(record.first, (template_segment + 1) % 2, coverage_block,
                                                         target)) {
                                return false;
                            }
                            if (!EvalReferenceStatistics(record.second, template_segment, coverage_block, target)) {
                                return false;
                            }
                        }
                    }
                }
            }
        }
    }

    return true;
}

bool BamIngestionEngine::SignsOfPairsWithNamesNotIdentical(const DataStats& target) {
    // All reads have been processed, so first_read_records_ must be empty
    if (first_read_records_.size()) {
        printErr
            << "Did not find a match for all paired reads: Read names of reads from the same pair are not identical or "
               "reads from a pair are missing. Cannot continue.\nExample position (ReferenceSequence:StartPosition): "
            << target.reference_->ReferenceIdFirstPart((*first_read_records_.begin())->record_.rID) << ":"
            << (*first_read_records_.begin())->record_.beginPos
            << "\nExample read name: " << (*first_read_records_.begin())->record_.qName << std::endl;
        return true;
    }

    return false;
}

void BamIngestionEngine::PrepareReadIn(uintQual size_mapping_quality, uintReadLen size_indel,
                                       uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads, DataStats& target) {
    uintReadLen size_pos = max(target.read_lengths_.at(0).to(), target.read_lengths_.at(1).to());

    uintFragCount num_reads(0);
    uintNucCount num_bases(0);
    for (int template_segment = kTemplateSegments; template_segment--;) {
        for (auto len = target.read_lengths_.at(template_segment).from();
             len < target.read_lengths_.at(template_segment).to(); ++len) {
            num_reads += target.read_lengths_.at(template_segment).at(len);
            num_bases += target.read_lengths_.at(template_segment).at(len) * len;
        }
    }

    target.adapters_.PrepareAdapters(size_pos, target.phred_quality_offset_);
    target.coverage_.Prepare(Divide(num_bases, target.reference_->TotalSize()), Divide(num_bases, num_reads),
                             target.maximum_read_length_on_reference_);
    target.duplicates_.PrepareTmpDuplicationVector();
    target.errors_.Prepare(target.tiles_.NumTiles(), target.maximum_quality_ + 1, size_pos, size_indel);
    target.fragment_distribution_.Prepare(*target.reference_, target.maximum_insert_length_, max_ref_seq_bin_size,
                                          reads_per_frag_len_bin_, lowq_reads_per_frag_len_bin_, num_threads);
    reads_per_frag_len_bin_.clear();
    reads_per_frag_len_bin_.shrink_to_fit();
    lowq_reads_per_frag_len_bin_.clear();
    lowq_reads_per_frag_len_bin_.shrink_to_fit();
    target.qualities_.Prepare(target.tiles_.NumTiles(), target.maximum_quality_ + 1, size_pos,
                              target.maximum_insert_length_);

    // Prepare vector in this class
    target.read_sequence_stats_.PrepareAccumulators(size_mapping_quality, size_pos);

    for (auto template_segment = kTemplateSegments; template_segment--;) {
        SetDimensions(target.tmp_read_lengths_by_fragment_length_.at(template_segment),
                      target.maximum_insert_length_ + 1, size_pos);
        SetDimensions(target.tmp_non_mapped_read_lengths_by_fragment_length_.at(template_segment),
                      target.maximum_insert_length_ + 1, size_pos);

        target.tmp_gc_read_content_reference_.at(template_segment).resize(kGCBins);
        target.tmp_gc_read_content_mapped_.at(template_segment).resize(kGCBins);

        for (auto strand = kStrands; strand--;) {
            for (auto base = kNumBases; base--;) {
                target.tmp_sequence_content_reference_.at(template_segment)
                    .at(strand)
                    .at(base)
                    .resize(target.maximum_read_length_on_reference_ + 1);
            }
        }
    }
}

void BamIngestionEngine::FinishReadIn(DataStats& target) {
    target.adapters_.Finalize(); // Collect ambigous adapters into a single one
    target.errors_.Finalize();
    target.fragment_distribution_.Finalize();
    target.qualities_.Finalize(SumVect(target.read_lengths_.at(0)));
    target.read_sequence_stats_.Finalize();

    // Copy vectors to final ones
    for (auto template_segment = kTemplateSegments; template_segment--;) {
        target.non_mapped_read_lengths_by_fragment_length_.at(template_segment)
            .Acquire(target.tmp_non_mapped_read_lengths_by_fragment_length_.at(template_segment));
        ShrinkVect(target.non_mapped_read_lengths_by_fragment_length_.at(template_segment));
        for (auto frag_len = target.non_mapped_read_lengths_by_fragment_length_.at(template_segment).from();
             frag_len < target.non_mapped_read_lengths_by_fragment_length_.at(template_segment).to(); ++frag_len) {
            for (auto read_len =
                     target.non_mapped_read_lengths_by_fragment_length_.at(template_segment).at(frag_len).from();
                 read_len < target.non_mapped_read_lengths_by_fragment_length_.at(template_segment).at(frag_len).to();
                 ++read_len) {
                target.tmp_read_lengths_by_fragment_length_.at(template_segment).at(frag_len).at(read_len) +=
                    target.non_mapped_read_lengths_by_fragment_length_.at(template_segment).at(frag_len).at(read_len);
            }
        }
        target.read_lengths_by_fragment_length_.at(template_segment)
            .Acquire(target.tmp_read_lengths_by_fragment_length_.at(template_segment));

        target.gc_read_content_reference_.at(template_segment)
            .Acquire(target.tmp_gc_read_content_reference_.at(template_segment));
        target.gc_read_content_mapped_.at(template_segment)
            .Acquire(target.tmp_gc_read_content_mapped_.at(template_segment));

        for (auto strand = kStrands; strand--;) {
            for (auto base = kNumBases; base--;) {
                target.sequence_content_reference_.at(template_segment)
                    .at(strand)
                    .at(base)
                    .Acquire(target.tmp_sequence_content_reference_.at(template_segment).at(strand).at(base));
            }
        }
    }
}

void BamIngestionEngine::Shrink(DataStats& target) {
    target.Shrink();
}

// Deactivation of bias calculation only for speeding up of tests
bool BamIngestionEngine::Calculate(uintNumThreads num_threads, DataStats& target) {
    // Calculate the duplicates from observation of fragments at same position and reference characteristics
    if (!target.fragment_distribution_.FinalizeBiasCalculation(*target.reference_, num_threads, target.duplicates_)) {
        return false;
    }

    // Calculate coverage corrected for low quality sites
    uintFragCount num_reads(0);
    uintNucCount num_bases(0);
    for (int template_segment = kTemplateSegments; template_segment--;) {
        for (auto len = target.read_lengths_.at(template_segment).from();
             len < target.read_lengths_.at(template_segment).to(); ++len) {
            num_reads += target.read_lengths_.at(template_segment).at(len);
            num_bases += target.read_lengths_.at(template_segment).at(len) * len;
        }
    }
    target.corrected_coverage_ = target.fragment_distribution_.CorrectedCoverage(
        *target.reference_, static_cast<double>(num_bases) / num_reads,
        target.coverage_.Coverage()); // Must be called after FinalizeBiasCalculation, where the coverage is corrected
    printInfo << "Estimated a corrected mean coverage of " << target.corrected_coverage_ << "x" << std::endl;

    return true;
}

bool BamIngestionEngine::OrderOfBamFileCorrect(const seqan::BamAlignmentRecord& record,
                                               pair<uintRefSeqId, uintSeqLen> last_record_pos) {
    if (!hasFlagUnmapped(record)) {
        if (record.rID < last_record_pos.first ||
            (record.rID == last_record_pos.first && record.beginPos < last_record_pos.second)) {
            printErr << "Bamfile is not sorted by position. Detected at position " << record.beginPos
                     << " with read name: " << record.qName << std::endl;
            return false;
        }
        last_record_pos = {record.rID, record.beginPos};
    }

    return true;
}

bool BamIngestionEngine::PreRun(BamFileIn& bam, const char* bam_file, BamHeader& header, uintQual& size_mapping_quality,
                                uintReadLen& size_indel, DataStats& target) {
    printInfo << "Starting PreRun" << std::endl;

    bool error = false;
    BamAlignmentRecord record;
    pair<uintRefSeqId, uintSeqLen> last_record_pos{0, 0};

    uintReadLen read_length_on_reference;
    uintQual min_mapq(255), max_mapq(0);
    size_indel = 0;

    uintRefLenCalc num_bins(0);
    for (auto ref_seq = target.reference_->NumberSequences(); ref_seq--;) {
        num_bins += target.reference_->SequenceLength(ref_seq) / target.maximum_insert_length_ + 1;
    }
    reads_per_frag_len_bin_.resize(num_bins);
    lowq_reads_per_frag_len_bin_.resize(num_bins);
    uintRefSeqId cur_ref_seq(0);
    uintRefLenCalc ref_seq_start_bin(0);

    target.adapters_.PrepareAdapterPrediction();
    do {
        try {
            ++target.total_number_reads_; // Counter total_number_reads_ for read in, later in Calculate correct it with
                                          // the number of accepted reads from read_length_

            readRecord(record, bam);

            if (target.IsValidRecord(record)) {
                if (OrderOfBamFileCorrect(record, last_record_pos)) {
                    if (!hasFlagSecondary(record) && !hasFlagSupplementary(record)) { // Ignore supplementary reads
                        if (hasFlagFirst(record)) {
                            target.tiles_.EnterTile(record.qName); // Names of records in a pair are identical, so tile
                                                                   // must only be identified once

                            ++target.read_lengths_.at(0)[length(record.qual)];
                        } else {
                            ++target.read_lengths_.at(1)[length(record.qual)];
                        }

                        // Determine read length range on reference
                        if (PotentiallyValid(record, target)) {
                            read_length_on_reference = GetReadLengthOnReference(record, size_indel);
                            SetToMax(target.maximum_read_length_on_reference_, read_length_on_reference);
                            SetToMin(target.minimum_read_length_on_reference_, read_length_on_reference);

                            while (record.rID > cur_ref_seq) {
                                ref_seq_start_bin +=
                                    target.reference_->SequenceLength(cur_ref_seq) / target.maximum_insert_length_ + 1;
                                ++cur_ref_seq;
                            }
                            ++reads_per_frag_len_bin_.at(ref_seq_start_bin +
                                                         record.beginPos / target.maximum_insert_length_);
                        }
                        if (!hasFlagUnmapped(record) && record.mapQ < target.minimum_mapping_quality_) {
                            while (record.rID > cur_ref_seq) {
                                ref_seq_start_bin +=
                                    target.reference_->SequenceLength(cur_ref_seq) / target.maximum_insert_length_ + 1;
                                ++cur_ref_seq;
                            }
                            ++lowq_reads_per_frag_len_bin_.at(ref_seq_start_bin +
                                                              record.beginPos / target.maximum_insert_length_);
                        }

                        // Determine quality range
                        for (auto i = length(record.qual); i--;) {
                            SetToMax(target.maximum_quality_, at(record.qual, i));
                            SetToMin(target.minimum_quality_, at(record.qual, i));
                        }

                        SetToMax(max_mapq, record.mapQ);
                        SetToMin(min_mapq, record.mapQ);

                        target.adapters_.ExtractAdapterPart(record);
                    }
                } else {
                    error = true;
                }
            } else {
                error = true;
            }
        } catch (const Exception& e) {
            error = true;
            printErr << "Could not read record " << target.total_number_reads_ << " in " << bam_file << ": " << e.what()
                     << std::endl;
        }
    } while (!atEnd(bam));

    if (error) {
        return false;
    }

    // Jump back to beginning of bam file
    seqan::close(bam);
    if (!open(bam, bam_file)) {
        printErr << "Could not open " << bam_file << " for reading." << std::endl;
        return false;
    }

    try {
        readHeader(header, bam);
    } catch (const Exception& e) {
        printErr << "Could not read header in " << bam_file << ": " << e.what() << std::endl;
        return false;
    }

    for (uintTempSeq template_segment = kTemplateSegments; template_segment--;) {
        target.read_lengths_.at(template_segment).Shrink();
    }

    if (target.minimum_quality_ < kPhredIlluminaOffset) {
        if (target.minimum_quality_ < kPhredSangerOffset) {
            printErr << "Minimum quality is " << target.minimum_quality_
                     << ", which is lower than the 33 from Sanger encoding. Are those values correct?" << std::endl;
            return false;
        }
        target.phred_quality_offset_ = kPhredSangerOffset; // Sanger format
    } else {
        target.phred_quality_offset_ = kPhredIlluminaOffset; // Old Illumina format
    }

    target.minimum_quality_ -= target.phred_quality_offset_;
    target.maximum_quality_ -= target.phred_quality_offset_;

    printInfo << "Finished PreRun\nTotal number of reads: " << target.total_number_reads_ << '\n'
              << "Read length: " << min(target.read_lengths_.at(0).from(), target.read_lengths_.at(1).from()) << " - "
              << (max(target.read_lengths_.at(0).to(), target.read_lengths_.at(1).to()) - 1) << '\n'
              << "Read length on reference: " << target.minimum_read_length_on_reference_ << " - "
              << target.maximum_read_length_on_reference_ << '\n'
              << "Quality: " << static_cast<uintQualPrint>(target.minimum_quality_) << " - "
              << static_cast<uintQualPrint>(target.maximum_quality_) << " ( "
              << static_cast<uintQualPrint>(target.phred_quality_offset_) << "-based encoding )\n"
              << "Mapping Quality: " << static_cast<uintQualPrint>(min_mapq) << " - "
              << static_cast<uintQualPrint>(max_mapq) << '\n'
              << "Maximum InDel Length: " << size_indel << std::endl;

    size_mapping_quality = max_mapq + 1;
    ++size_indel; // Add one to go from index of last value to size

    return target.adapters_.PredictAdapters();
}

bool BamIngestionEngine::ReadRecords(BamFileIn& bam, bool& not_done, ThreadData& thread_data, DataStats& target) {
    CoverageStats::FullRecord* record;
    CoverageStats::CoverageBlock* cov_block(nullptr);
    lock_guard<mutex> lock(read_mutex_);
    try {
        while (thread_data.rec_store_.size() < kBatchSize && !atEnd(bam) && reading_success_) {
            ++read_records_;
            record = new CoverageStats::FullRecord;
            readRecord(record->record_, bam);

            if (hasFlagSecondary(record->record_) || hasFlagSupplementary(record->record_)) {
                delete record;
            } else {
                if (PotentiallyValid(record->record_, target)) {
                    // Use record->record_.beginPos-maximum_read_length_on_reference_ as start, because of potentially
                    // soft-clipped bases at the beginning
                    if (!target.coverage_.EnsureSpace(
                            record->record_.rID,
                            (record->record_.beginPos > target.maximum_read_length_on_reference_
                                 ? record->record_.beginPos - target.maximum_read_length_on_reference_
                                 : 0),
                            record->record_.beginPos + target.maximum_read_length_on_reference_, record,
                            *target.reference_)) {
                        return false;
                    }
                }

                // Add low q sites
                if (!target.QualitySufficient(record->record_) && !hasFlagUnmapped(record->record_) &&
                    !target.reference_->ReferenceSequenceExcluded(record->record_.rID) && // low quality and not
                                                                                          // excluded
                    (hasFlagNextUnmapped(record->record_) || record->record_.rID != record->record_.rNextId ||
                     (!hasFlagRC(record->record_) &&
                      (record->record_.beginPos < record->record_.pNext ||
                       record->record_.beginPos > record->record_.pNext + length(record->record_.seq))) ||
                     (hasFlagRC(record->record_) && (record->record_.beginPos > record->record_.pNext ||
                                                     record->record_.beginPos + length(record->record_.seq) <
                                                         record->record_.pNext)))) { // not in an adapter pair
                    uintSeqLen start_pos, end_pos;
                    target.GetReadPosOnReference(start_pos, end_pos, record->record_);

                    if (!target.reference_->FragmentExcluded(thread_data.last_exclusion_region_id_,
                                                             thread_data.last_exclusion_ref_seq_, record->record_.rID,
                                                             start_pos, end_pos)) {
                        if (hasFlagRC(record->record_)) {
                            target.fragment_distribution_.AddLowQSiteEnd(
                                record->record_.rID, end_pos, *target.reference_, thread_data.fragment_distribution_);
                        } else {
                            target.fragment_distribution_.AddLowQSiteStart(
                                record->record_.rID, start_pos, *target.reference_, thread_data.fragment_distribution_);
                        }
                    }
                }

                // Work with the pairs and not individual reads, so first combine reads to pairs, by storing first and
                // then processing both reads at the same time
                CoverageStats::FullRecord* record_first;
                if (IsSecondRead(record, record_first, cov_block, target)) {
                    if (length(record_first->record_.qual) == length(record_first->record_.seq) &&
                        length(record->record_.qual) == length(record->record_.seq)) {
                        thread_data.rec_store_.emplace_back(record_first, record);
                    } else {
                        printWarn << "Ignoring read pair, because sequence and quality have different length for "
                                  << record->record_.qName << ": Seq(" << length(record_first->record_.seq) << ", "
                                  << length(record->record_.seq) << ") Qual(" << length(record_first->record_.qual)
                                  << ", " << length(record->record_.qual) << ')' << std::endl;
                    }
                }
            }

            if (target.total_number_reads_ > kBatchSize && !(read_records_ % (target.total_number_reads_ / 20))) {
                lock_guard<mutex> lock(print_mutex_);
                printInfo << "Read "
                          << static_cast<uintPercentPrint>(Percent(read_records_, target.total_number_reads_))
                          << "% of the reads." << std::endl;
            }
        }
        not_done = !atEnd(bam);
    } catch (const Exception& e) {
        lock_guard<mutex> lock(print_mutex_);
        printErr << "Could not read record " << read_records_ << ": " << e.what() << std::endl;
        return false;
    }

    return true;
}

void BamIngestionEngine::ReadThread(BamIngestionEngine& engine, BamFileIn& bam, size_t thread_idx, DataStats& target) {
    CoverageStats::CoverageBlock* cov_block;
    uintFragCount processed_fragments(0);
    bool not_done(true);
    auto first_ref_seq = target.reference_->FirstNotExcludedSequence();
    if (first_ref_seq >= target.reference_->NumberSequences()) {
        printErr << "No contig longer than " << target.maximum_insert_length_ << std::endl;
        engine.reading_success_ = false;
        first_ref_seq = 0; // Just so that the next command does not crash, before we shut down
    }
    ThreadData thread_data(target.maximum_insert_length_, target.reference_->StartExclusion(first_ref_seq));
    thread_data.rec_store_.reserve(engine.kBatchSize);
    uintRefSeqId still_needed_reference_sequence(0);
    uintSeqLen still_needed_position(0);
    while (not_done && engine.reading_success_) {
        // ReadRecords:
        // Reads in records and bundles them to batches of paired records
        // Registers reads as unprocessed in first coverage block, to make sure blocks are not removed before all reads
        // in them are handled Adds pointers of potentially valid reads to the last block they potentially overlap (read
        // length on reference yet unknown), so the coverage part can be handled after all coverage information are
        // gathered
        if (engine.ReadRecords(bam, not_done, thread_data, target)) {
            cov_block = nullptr;
            // Process batch
            for (auto& rec : thread_data.rec_store_) {
                // EvalRecord:
                // Apply multiple filters to check whether fragments are valid
                // Set to_ref_pos for valid reads, that to_ref_pos==0 is used in coverage_.CleanUp to remove the
                // non-valid reads Calculate all the statistics that are independent of other reads Fills the coverage
                // information in the coverage blocks
                if (engine.EvalRecord(rec, thread_data, target)) {
                    if (engine.PotentiallyValidFirst(rec.first->record_, target)) {
                        // coverage_.RemoveFragment
                        // Marks the reads as processed in the first block of first read, so the fully processed
                        // coverage blocks can be handled after all reads are processed Use
                        // record->record_.beginPos-maximum_read_length_on_reference_ as start, because of potentially
                        // soft-clipped bases at the beginning
                        target.coverage_.RemoveFragment(
                            rec.first->record_.rID,
                            (rec.first->record_.beginPos > target.maximum_read_length_on_reference_
                                 ? rec.first->record_.beginPos - target.maximum_read_length_on_reference_
                                 : 0),
                            cov_block, processed_fragments);
                    } else {
                        // Remove reads that have not been added to a coverage block
                        delete rec.first;
                        delete rec.second;
                    }
                } else {
                    engine.reading_success_ = false;
                }
            }

            if (cov_block) { // To make sure we don't have a whole batch of non-valid reads as coverage_.CleanUp already
                             // handled blocks at the beginning and we need at least one block still existing to keep
                             // track of the position
                // coverage_.CleanUp:
                // Once all reads in it have been handled processes the coverage blocks: The coverage itself and all the
                // read statistics that depend on the other reads (systematic error) Removes the processed blocks
                still_needed_reference_sequence = target.coverage_.CleanUp(
                    still_needed_position, *target.reference_, target.qualities_, target.errors_,
                    target.phred_quality_offset_, cov_block, processed_fragments, thread_data.coverage_);
            }

            target.coverage_.PreLoadVariants(*target.reference_);

            target.fragment_distribution_.HandleReferenceSequencesUntil(
                still_needed_reference_sequence, still_needed_position, thread_data.fragment_distribution_,
                *target.reference_, target.duplicates_, engine.print_mutex_, thread_idx);
        } else {
            engine.reading_success_ = false;
        }

        thread_data.rec_store_.clear();
    }

    if (engine.reading_success_) {
        if (0 == --engine.running_threads_) {
            {
                lock_guard<mutex> lock(engine.finish_threads_mutex_);
                engine.finish_threads_ = true;
            }
            engine.finish_threads_cv_.notify_all();

            // This finalization might take a bit of time so do it already here
            if (!target.coverage_.Finalize(*target.reference_, target.qualities_, target.errors_,
                                           target.phred_quality_offset_, engine.print_mutex_, thread_data.coverage_)) {
                engine.reading_success_ = false;
            }
        } else {
            unique_lock<mutex> lock(engine.finish_threads_mutex_);
            engine.finish_threads_cv_.wait(lock, [&engine] { return engine.finish_threads_; });
        }

        if (engine.reading_success_) {
            target.fragment_distribution_.FinishThreads(thread_data.fragment_distribution_, *target.reference_,
                                                        target.duplicates_, engine.print_mutex_, thread_idx);
            target.fragment_distribution_.AddThreadData(thread_data.fragment_distribution_);
        }
    }
}

reseq::uintReadLen BamIngestionEngine::GetReadLengthOnReference(const BamAlignmentRecord& record) {
    uintReadLen real_read_length = length(record.seq);

    // Update real_read_length based on insertion and deletion
    for (const auto& cigar_element : record.cigar) {
        switch (cigar_element.operation) {
        case 'N':
        case 'D':
            real_read_length += cigar_element.count;
            break;
        case 'I':
            real_read_length -= cigar_element.count;
            break;
        default:
            break;
        }
    }
    return real_read_length;
}

reseq::uintReadLen BamIngestionEngine::GetReadLengthOnReference(const BamAlignmentRecord& record,
                                                                uintReadLen& max_indel) {
    uintReadLen real_read_length = length(record.seq);

    // Update real_read_length based on insertion and deletion
    for (const auto& cigar_element : record.cigar) {
        switch (cigar_element.operation) {
        case 'N':
        case 'D':
            SetToMax(max_indel, cigar_element.count);
            real_read_length += cigar_element.count;
            break;
        case 'I':
            SetToMax(max_indel, cigar_element.count);
            real_read_length -= cigar_element.count;
            break;
        default:
            break;
        }
    }
    return real_read_length;
}

bool BamIngestionEngine::Run(const char* bam_file, const char* adapter_file, const char* adapter_matrix,
                             const string& variant_file, uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads,
                             bool calculate_bias, DataStats& target) {
    // Initialize engine state
    reading_success_ = true;
    read_records_ = 0;
    reads_in_unmapped_pairs_without_adapters_ = 0;
    reads_in_unmapped_pairs_with_adapters_ = 0;
    reads_with_low_quality_with_adapters_ = 0;
    reads_with_low_quality_without_adapters_ = 0;
    reads_on_too_short_fragments_ = 0;
    reads_in_excluded_regions_ = 0;
    reads_used_ = 0;

    if (calculate_bias) {
        target.fragment_distribution_.ActivateBiasCalculation();
    } else {
        target.fragment_distribution_.DeactivateBiasCalculation(); // Only for speeding up tests
    }

    bool success = true;
    BamFileIn bam;

    if (target.reference_) {
        if (!open(bam, bam_file)) {
            printErr << "Could not open " << bam_file << " for reading." << std::endl;
            success = false;
        } else if (atEnd(bam)) {
            printErr << bam_file << " does not contain any sequences." << std::endl;
            success = false;
        } else {
            BamHeader header;
            try {
                readHeader(header, bam);
            } catch (const Exception& e) {
                printErr << "Could not read header in " << bam_file << ": " << e.what() << std::endl;
                success = false;
            }

            if (success) {
                const auto& bam_context = context(bam);

                if (length(contigNames(bam_context)) == target.reference_->NumberSequences()) {
                    bool reference_ordering_identical = true;
                    for (auto i = length(contigNames(bam_context)); i--;) {
                        if (!(target.reference_->ReferenceIdFirstPart(i) == at(contigNames(bam_context), i))) {
                            reference_ordering_identical = false;
                            printErr << "The ordering of the reference sequences in the bam file is not identical to "
                                        "the reference file\n"
                                     << at(contigNames(bam_context), i) << '\n'
                                     << target.reference_->ReferenceId(i) << std::endl;
                            success = false;
                            break;
                        }
                    }

                    if (reference_ordering_identical) {
                        uintQual size_mapping_quality;
                        uintReadLen size_indel;

                        target.reference_->PrepareExclusionRegions();
                        target.reference_->ObtainExclusionRegions(
                            target.reference_->NumberSequences(),
                            target.maximum_insert_length_); // At the end we need all of them together anyways, so it
                                                            // makes no sense to read them in one after another

                        if ((0 == strlen(adapter_file) && 0 == strlen(adapter_matrix)) ||
                            target.adapters_.LoadAdapters(adapter_file, adapter_matrix)) {
                            if (PreRun(bam, bam_file, header, size_mapping_quality, size_indel, target)) {
                                if (!variant_file.empty()) {
                                    if (target.reference_->PrepareVariantFile(variant_file)) {
                                        if (!target.reference_->ReadFirstVariantPositions()) {
                                            success = false;
                                        }
                                    } else {
                                        success = false;
                                    }
                                }

                                if (success) {
                                    PrepareReadIn(size_mapping_quality, size_indel, max_ref_seq_bin_size, num_threads,
                                                  target);

                                    printInfo << "Starting main read-in" << std::endl;

                                    running_threads_ = num_threads;
                                    finish_threads_ = false;
                                    {
                                        std::vector<std::jthread> threads;
                                        threads.reserve(num_threads);
                                        for (decltype(num_threads) i = 0; i < num_threads; ++i) {
                                            threads.emplace_back([this, &bam, i, &target](const std::stop_token&) {
                                                ReadThread(*this, bam, i, target);
                                            });
                                        }
                                        // jthread destructors join on scope exit
                                    }

                                    if (reading_success_) {
                                        printInfo << "Finished processing all reads." << std::endl;
                                    } else {
                                        success = false;
                                    }
                                }
                            } else {
                                success = false;
                            }
                        } else {
                            success = false;
                        }
                    }
                } else {
                    printErr << "The number of the reference sequences in the bam file is not identical to the number "
                                "of sequences in the reference file"
                             << std::endl;
                    success = false;
                }
            }
        }
    } else {
        printErr << "No reference provided" << std::endl;
        success = false;
    }

    seqan::close(bam);

    target.reference_->ClearAllVariantPositions();

    if (!success) {
        target.total_number_reads_ = 0; // Marking that the reading in was not successful
        return false;
    }

    uintFragCount reads_in_wrong_place = target.total_number_reads_ - reads_in_unmapped_pairs_without_adapters_ -
                                         reads_in_unmapped_pairs_with_adapters_ -
                                         reads_with_low_quality_with_adapters_ -
                                         reads_with_low_quality_without_adapters_ - reads_on_too_short_fragments_ -
                                         reads_in_excluded_regions_ - reads_used_;
    printInfo << "Of the " << target.total_number_reads_ << " reads in the file" << std::endl;
    printInfo << reads_used_ << " (" << static_cast<uintPercentPrint>(Percent(reads_used_, target.total_number_reads_))
              << "\%) could be used for all statistics" << std::endl;
    printInfo << reads_with_low_quality_with_adapters_ << " ("
              << static_cast<uintPercentPrint>(
                     Percent(reads_with_low_quality_with_adapters_, target.total_number_reads_))
              << "\%) had too low mapping qualities with adapters detected" << std::endl;
    printInfo << reads_with_low_quality_without_adapters_ << " ("
              << static_cast<uintPercentPrint>(
                     Percent(reads_with_low_quality_without_adapters_, target.total_number_reads_))
              << "\%) had too low mapping qualities without adapters detected" << std::endl;
    printInfo << reads_in_excluded_regions_ << " ("
              << static_cast<uintPercentPrint>(Percent(reads_in_excluded_regions_, target.total_number_reads_))
              << "\%) mapped to excluded regions of the reference" << std::endl;
    printInfo << reads_in_wrong_place << " ("
              << static_cast<uintPercentPrint>(Percent(reads_in_wrong_place, target.total_number_reads_))
              << "\%) were mapping too far apart from their partner or in wrong direction" << std::endl;
    printInfo << reads_on_too_short_fragments_ << " ("
              << static_cast<uintPercentPrint>(Percent(reads_on_too_short_fragments_, target.total_number_reads_))
              << "\%) are on reference sequences that are too short" << std::endl;
    printInfo << reads_in_unmapped_pairs_with_adapters_ << " ("
              << static_cast<uintPercentPrint>(
                     Percent(reads_in_unmapped_pairs_with_adapters_, target.total_number_reads_))
              << "\%) were in an unmapped pair with adapters detected" << std::endl;
    printInfo << reads_in_unmapped_pairs_without_adapters_ << " ("
              << static_cast<uintPercentPrint>(
                     Percent(reads_in_unmapped_pairs_without_adapters_, target.total_number_reads_))
              << "\%) were in an unmapped pair without adapters detected" << std::endl;

    if (0 == reads_used_ &&
        calculate_bias) { // In case we are testing and not calculating the bias afterwards, we can continue running
        printErr << "No reads in the file passed all criteria to be used for the statistics" << std::endl;
        success = false;
    }

    if (!success || SignsOfPairsWithNamesNotIdentical(target)) {
        target.total_number_reads_ = 0; // Marking that the reading in was not successful
        return false;
    }
    FinishReadIn(target);

    Shrink(target);
    if (!Calculate(num_threads, target)) {
        target.total_number_reads_ = 0; // Marking that the reading in was not successful
        return false;
    }

    target.reference_
        ->ClearAllExclusionRegions(); // Don't remove any of it earlier because we still need all of them in Calculate

    target.creation_time_ =
        std::chrono::time_point_cast<std::chrono::seconds>(std::chrono::system_clock::now()).time_since_epoch().count();

    return true;
}
