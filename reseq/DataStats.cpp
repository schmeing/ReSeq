#include "DataStats.h"
using reseq::DataStats;
using reseq::Vect;

#include <exception>
using std::exception;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <limits>
using std::numeric_limits;
// include <string>
using std::string;

#include "archive_format.h"
#include "BamIngestionEngine.h"
#include "logging.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

// include <seqan/bam_io.h>
using seqan::BamAlignmentRecord;

// include "utilities.hpp"
using reseq::utilities::at;
using reseq::utilities::CreateDir;
using reseq::utilities::FileExists;

DataStats::DataStats(Reference* ref, uintSeqLen maximum_insert_length, uintQual minimum_mapping_quality)
    : reference_(ref), maximum_insert_length_(maximum_insert_length),
      minimum_mapping_quality_(
          minimum_mapping_quality), // For arabidopsis thaliana mapping with bowtie2 a strange clustering of fragments
                                    // at same positions was observed for lower qualities than 10
      minimum_quality_(255), maximum_quality_(0), minimum_read_length_on_reference_(numeric_limits<uintReadLen>::max()),
      maximum_read_length_on_reference_(0), total_number_reads_(0) {
    for (uintTempSeq template_segment = kTemplateSegments; template_segment--;) {
        read_lengths_.at(template_segment).SetOffset(1); // As a read with length of 0 cannot be considered a read, this
                                                         // length can be excluded right from the start
    }
}

bool DataStats::IsValidRecord(const BamAlignmentRecord& record) {
    if (!hasFlagMultiple(record)) {
        printErr << "Read '" << record.qName << "' is not paired." << std::endl;
        return false;
    }
    if (hasFlagFirst(record) == hasFlagLast(record)) {
        printErr << "Read '" << record.qName << "' is either first and second in its template or none of it."
                 << std::endl;
        return false;
    }
    if (0 == length(record.seq)) {
        printErr << "Read '" << record.qName << "' does not have the sequence stored." << std::endl;
        return false;
    }
    if (length(record.seq) != length(record.qual)) {
        printErr << "Read '" << record.qName << "' does not have the same length for the sequence and its quality."
                 << std::endl;
        return false;
    }

    return true;
}

reseq::uintReadLen DataStats::GetReadLengthOnReference(const BamAlignmentRecord& record) {
    return BamIngestionEngine::GetReadLengthOnReference(record);
}

reseq::uintReadLen DataStats::GetReadLengthOnReference(const BamAlignmentRecord& record, uintReadLen& max_indel) {
    return BamIngestionEngine::GetReadLengthOnReference(record, max_indel);
}

void DataStats::GetReadPosOnReference(uintSeqLen& start_pos, uintSeqLen& end_pos,
                                      const BamAlignmentRecord& record) const {
    end_pos = GetReadLengthOnReference(record); // start_pos will be added after it is corrected for soft-clipping
    start_pos = record.beginPos;

    if ('S' == at(record.cigar, 0).operation) {
        if (at(record.cigar, 0).count < start_pos) {
            start_pos -= at(record.cigar, 0).count;
        } else {
            start_pos = 0;
        }
    } else if ('H' == at(record.cigar, 0).operation) {
        end_pos -= at(record.cigar, 0).count;

        if (2 <= length(record.cigar) && 'S' == at(record.cigar, 1).operation) {
            if (at(record.cigar, 1).count < start_pos) {
                start_pos -= at(record.cigar, 1).count;
            } else {
                start_pos = 0;
            }
        }
    }

    end_pos += start_pos;

    if ('H' == at(record.cigar, length(record.cigar) - 1).operation && 2 <= length(record.cigar)) {
        end_pos -= at(record.cigar, length(record.cigar) - 1).count;
    }

    if (reference_->SequenceLength(record.rID) < end_pos) {
        end_pos = reference_->SequenceLength(record.rID); // Can occur due to soft-clipping at the end
    }

    return;
}

bool DataStats::ReadBam(const char* bam_file, const char* adapter_file, const char* adapter_matrix,
                        const string& variant_file, uintSeqLen max_ref_seq_bin_size, uintNumThreads num_threads,
                        bool calculate_bias) {
    BamIngestionEngine engine;
    return engine.Run(bam_file, adapter_file, adapter_matrix, variant_file, max_ref_seq_bin_size, num_threads,
                      calculate_bias, *this);
}

bool DataStats::Load(const char* archive_file) {
    if (!FileExists(archive_file)) {
        printErr << "File '" << archive_file << "' does not exists or no read permission given." << std::endl;
        return false;
    }

    try {
        ifstream ifs(archive_file, std::ios::binary);
        auto fmt = reseq::format::DetectFormat(ifs, reseq::format::kStatsMagic);

        switch (fmt) {
        case reseq::format::ArchiveFormat::kCompressedBinaryV1: {
            boost::iostreams::filtering_istream fis;
            fis.push(boost::iostreams::gzip_decompressor());
            fis.push(ifs);
            boost::archive::binary_iarchive ia(fis);
            ia >> *this;
            break;
        }
        case reseq::format::ArchiveFormat::kText: {
            boost::archive::text_iarchive ia(ifs);
            ia >> *this;
            break;
        }
        case reseq::format::ArchiveFormat::kUnsupportedBinaryVersion:
            printErr << "Unsupported profile format version in '" << archive_file << "'. Please upgrade ReSeq."
                     << std::endl;
            return false;
        }
    } catch (const exception& e) {
        printErr << "Could not load data statistics: " << e.what() << std::endl;
        return false;
    }

    return true;
}

bool DataStats::Save(const char* archive_file, bool text_format) const {
    try {
        CreateDir(archive_file);

        ofstream ofs(archive_file, std::ios::binary);
        if (text_format) {
            boost::archive::text_oarchive oa(ofs);
            oa << *this;
        } else {
            reseq::format::WriteHeader(ofs, reseq::format::kStatsMagic);
            boost::iostreams::filtering_ostream fos;
            fos.push(boost::iostreams::gzip_compressor());
            fos.push(ofs);
            boost::archive::binary_oarchive oa(fos);
            oa << *this;
            fos.flush();
        }
    } catch (const exception& e) {
        printErr << "Could not save data statistics: " << e.what() << std::endl;
        return false;
    }

    return true;
}

void DataStats::Shrink() {
    adapters_.Shrink();
    coverage_.Shrink();
    errors_.Shrink();
    qualities_.Shrink();
    tiles_.Shrink();

    for (uintTempSeq template_segment = kTemplateSegments; template_segment--;) {
        ShrinkVect(read_lengths_by_fragment_length_.at(template_segment));
    }
}

void DataStats::PrepareGeneral() {
    // Set total_number_reads_ to the number of accepted reads as it is the number of read records after read in
    total_number_reads_ = SumVect(read_lengths_.at(0)) + SumVect(read_lengths_.at(1));

    adapters_.SumCounts();
}

void DataStats::PrepareProcessing() {
    PrepareGeneral();

    adapters_.PrepareSimulation();
    qualities_.PrepareEstimation();
    errors_.PrepareSimulation();
}

void DataStats::PreparePlotting() {
    PrepareGeneral();
    coverage_.PreparePlotting();
    fragment_distribution_.PreparePlotting();
    errors_.PreparePlotting();
    qualities_.PreparePlotting();
}

void DataStats::PrepareTesting() {
    PrepareGeneral();

    coverage_.PreparePlotting();
    errors_.PreparePlotting();
    qualities_.PrepareTesting();
}
