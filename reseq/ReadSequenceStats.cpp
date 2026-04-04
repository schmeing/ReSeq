#include "ReadSequenceStats.h"
using reseq::ReadSequenceStats;

#include <array>
using std::array;

#include <seqan/bam_io.h>
using seqan::BamAlignmentRecord;
using seqan::Dna5;

#include "SeqQualityStats.hpp"
// include "utilities.hpp"
using reseq::utilities::at;
using reseq::utilities::Complement;
using reseq::utilities::Percent;

void ReadSequenceStats::PrepareAccumulators(uintQual size_mapping_quality, uintReadLen size_pos) {
    tmp_proper_pair_mapping_quality_.resize(size_mapping_quality);
    tmp_improper_pair_mapping_quality_.resize(size_mapping_quality);
    tmp_single_read_mapping_quality_.resize(size_mapping_quality);

    for (auto template_segment = kTemplateSegments; template_segment--;) {
        tmp_gc_read_content_.at(template_segment).resize(kGCBins);
        tmp_n_content_.at(template_segment).resize(101);
        for (auto base = kNumBasesN; base--;) {
            tmp_sequence_content_.at(template_segment).at(base).resize(size_pos);
        }
    }

    for (auto base = kNumBasesN; base--;) {
        tmp_homopolymer_distribution_.at(base).resize(size_pos);
    }
}

void ReadSequenceStats::EvalBaseLevelStats(CoverageStats::FullRecord* full_record, uintTempSeq template_segment,
                                           uintTempSeq strand, uintTileId tile_id, uintQual& paired_seq_qual,
                                           QualityStats& qualities, uintQual phred_quality_offset,
                                           uintQual maximum_quality) {
    const BamAlignmentRecord& record(full_record->record_);

    uintSeqLen pos_reversed(length(record.seq));

    SeqQualityStats<uintNucCount> seq_qual_stats;
    seq_qual_stats[maximum_quality]; // Resize to maximum quality
    for (auto qual : record.qual) {
        ++seq_qual_stats.at(qual - phred_quality_offset);
    }
    seq_qual_stats.Calculate();
    seq_qual_stats.CalculateProbabilityMean();
    full_record->sequence_quality_ = seq_qual_stats.mean_;

    Dna5 base;
    uintQual qual, last_qual(1);
    array<uintSeqLen, 5> read_bases = {0, 0, 0, 0, 0};
    uintReadLen homoquality_length = 0;
    uintReadLen homopolymer_length = 0;
    Dna5 homopolymer_nucleotide = 0;

    for (uintReadLen pos = 0; pos < length(record.seq); ++pos) {
        // In case of reversed sequences base and qual are read in reverse from bamfile so pos represents the position
        // in the fq file
        if (hasFlagRC(record)) {
            base = Complement::Dna5(at(record.seq, --pos_reversed));
            qual = at(record.qual, pos_reversed) - phred_quality_offset;
        } else {
            base = at(record.seq, pos);
            qual = at(record.qual, pos) - phred_quality_offset;
        }

        qualities.AddRawBase(template_segment, base, tile_id, strand, qual, seq_qual_stats.mean_, last_qual, pos);

        ++read_bases.at(base);
        ++tmp_sequence_content_.at(template_segment).at(base).at(pos);

        if (last_qual == qual) {
            ++homoquality_length;
        } else {
            if (pos) {
                qualities.AddRawHomoqualimer(last_qual, homoquality_length);
            }
            homoquality_length = 1;
        }

        if (homopolymer_nucleotide == base) {
            ++homopolymer_length;
        } else {
            ++tmp_homopolymer_distribution_.at(homopolymer_nucleotide).at(homopolymer_length);
            homopolymer_nucleotide = base;
            homopolymer_length = 1;
        }

        // Quality based on preceding quality
        last_qual = qual;
    }

    qualities.AddRawHomoqualimer(last_qual, homoquality_length);
    ++tmp_homopolymer_distribution_.at(homopolymer_nucleotide)
          .at(homopolymer_length); // Add the homopolymer at read end

    // Read level summaries of base level stats
    qualities.AddRawRead(paired_seq_qual, seq_qual_stats, template_segment, tile_id, read_bases, length(record.seq));

    ++tmp_gc_read_content_.at(template_segment)
          .at(Percent(read_bases.at(1) + read_bases.at(2), static_cast<uintSeqLen>(length(record.seq))));
    ++tmp_n_content_.at(template_segment).at(Percent(read_bases.at(4), static_cast<uintSeqLen>(length(record.seq))));
}

void ReadSequenceStats::IncrementMappingQuality(uintTempSeq quality_type, uintQual value) {
    // quality_type: 0 = proper pair, 1 = improper pair, 2 = single read
    switch (quality_type) {
    case 0:
        ++tmp_proper_pair_mapping_quality_.at(value);
        break;
    case 1:
        ++tmp_improper_pair_mapping_quality_.at(value);
        break;
    case 2:
        ++tmp_single_read_mapping_quality_.at(value);
        break;
    default:
        break;
    }
}

void ReadSequenceStats::Finalize() {
    proper_pair_mapping_quality_.Acquire(tmp_proper_pair_mapping_quality_);
    improper_pair_mapping_quality_.Acquire(tmp_improper_pair_mapping_quality_);
    single_read_mapping_quality_.Acquire(tmp_single_read_mapping_quality_);

    for (auto template_segment = kTemplateSegments; template_segment--;) {
        gc_read_content_.at(template_segment).Acquire(tmp_gc_read_content_.at(template_segment));
        n_content_.at(template_segment).Acquire(tmp_n_content_.at(template_segment));
        for (auto base = kNumBasesN; base--;) {
            sequence_content_.at(template_segment)
                .at(base)
                .Acquire(tmp_sequence_content_.at(template_segment).at(base));
        }
    }

    tmp_homopolymer_distribution_.at(0).at(0) = 0; // Remove homopolymers introduced by initialization values
    for (auto base = kNumBasesN; base--;) {
        homopolymer_distribution_.at(base).Acquire(tmp_homopolymer_distribution_.at(base));
    }
}

void ReadSequenceStats::Shrink() {
    // No additional shrinking needed for these histograms beyond what Acquire already does.
    // This method exists for consistency with the extraction pattern.
}
