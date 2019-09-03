#ifndef DATASTATSINTERFACE_H
#define DATASTATSINTERFACE_H

#include <map>
#include <stdint.h>
#include <string>
#include <utility>
#include <vector>

#include "DataStats.h"

namespace reseq{
	class DataStatsInterface{
	private:
		DataStats stats_;

	public:
		DataStatsInterface(Reference *ref);

		// Getter functions with std return values (used for python plotting)
		const char *AdapterName(uint16_t template_segment, uint16_t id) const;
		const std::vector<uint64_t> &AdapterCount(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &AdapterPolyATailLength() const;
		uint64_t AdapterOverrunBases(uint16_t nucleotide) const;

		uint32_t ErrorRatesByDistanceStart() const;
		uint32_t ErrorRatesByDistanceEnd() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ErrorRatesByDistance(uint32_t distance) const;
		uint16_t ErrorRatesByGCStart() const;
		uint16_t ErrorRatesByGCEnd() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ErrorRatesByGC(uint16_t gc) const;

		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &Coverage() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &CoverageStranded( bool reverse_strand ) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &CoverageStrandedPercent( bool reverse_strand ) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &CoverageStrandedPercentMinCov10( bool reverse_strand ) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &CoverageStrandedPercentMinCov20( bool reverse_strand ) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ErrorCoverage() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ErrorCoveragePercent() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ErrorCoveragePercentMinCov10() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ErrorCoveragePercentMinCov20() const;
		uint64_t ErrorCoveragePercentStranded( unsigned char forward_errors_percent, unsigned char reverse_errors_percent ) const;
		uint64_t ErrorCoveragePercentStrandedMinCov10( unsigned char forward_errors_percent, unsigned char reverse_errors_percent ) const;
		uint64_t ErrorCoveragePercentStrandedMinCov20( unsigned char forward_errors_percent, unsigned char reverse_errors_percent ) const;

		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ErrorsPerRead(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &CalledBasesByBaseQualityPerPreviousCalledBase(uint16_t template_segment, uint16_t ref_base, uint16_t called_base, uint16_t previous_base) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &CalledBasesByBaseQuality(uint16_t template_segment, uint16_t ref_base, uint16_t called_base) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &CalledBasesByPosition(uint16_t template_segment,uint16_t ref_base, uint16_t called_base) const;

		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &InDelErrorByLength(uint16_t indel_type) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &InDelErrorByPosition(uint16_t indel_type) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &InDelErrorByGC(uint16_t indel_type) const;

		const std::vector<double> &DispersionList() const;
		const std::vector<double> &MeanList() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &FragmentDuplicationNumber() const;

		const std::vector<uint64_t> &Abundance() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &InsertLengths() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &GCFragmentContent() const;
		const std::pair< std::vector<double>::size_type, std::vector<double> > &GCFragmentContentBias() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &OutskirtContent(uint16_t direction, uint16_t nucleotide) const;
		const std::vector<double> &FragmentSurroundingBiasByBase(uint16_t nucleotide) const;

		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityProbabilityMean(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityMinimum(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityFirstQuartile(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityMedian(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityThirdQuartile(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityMaximum(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityContent(uint16_t template_segment, unsigned char quality) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &HomoqualityDistribution(uint16_t quality) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &NucleotideQuality(uint16_t template_segment, uint16_t nucleotide) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &BaseQualityStatsReference(uint16_t template_segment, uint32_t read_position) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMeanReference(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMinimumReference(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityFirstQuartileReference(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMedianReference(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityThirdQuartileReference(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMaximumReference(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &AverageSequenceQualityForGC(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &BaseQualityStats(uint16_t template_segment, uint32_t read_position) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMean(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMinimum(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityFirstQuartile(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMedian(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityThirdQuartile(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMaximum(uint16_t template_segment) const;
		const std::pair< std::vector<signed char>::size_type, std::vector<signed char> > &TileQualityMeanDifference(uint16_t template_segment, uint16_t tile_id) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &BaseQualityMeanPerStrand(uint16_t strand) const;
		uint16_t BaseQualityForSequenceStart(uint16_t template_segment) const;
		uint16_t BaseQualityForSequenceEnd(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &BaseQualityForSequence(uint16_t template_segment, uint16_t seq_quality) const;
		uint16_t BaseQualityForPrecedingQualityStart(uint16_t template_segment) const;
		uint16_t BaseQualityForPrecedingQualityEnd(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &BaseQualityForPrecedingQuality(uint16_t template_segment, uint16_t seq_quality) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityMean(uint16_t template_segment) const;
		uint16_t SequenceQualityPairsStart() const;
		uint16_t SequenceQualityPairsEnd() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceQualityPairs(uint16_t seq_quality_first) const;
		const std::pair< std::vector<uint16_t>::size_type, std::vector<uint16_t> > &MeanSequenceQualityMeanByFragmentLength(uint16_t template_segment) const;
		const std::pair< std::vector<unsigned char>::size_type, std::vector<unsigned char> > &AverageSequenceQualityForBase(uint16_t template_segment, uint16_t nucleotide) const;

		const std::vector<uint16_t> &TileNames() const;
		const std::vector<uint64_t> &TileAbundance() const;

		unsigned char PhredQualityOffset() const;
		uint64_t TotalNumberReads() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ReadLengths(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ProperPairMappingQuality() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &ImproperPairMappingQuality() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SingleReadMappingQuality() const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &GCReadContent(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &GCReadContentReference(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &GCReadContentMapped(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &NContent(uint16_t template_segment) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceContent(uint16_t template_segment, uint16_t nucleotide) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &SequenceContentReference(uint16_t template_segment, uint16_t strand, uint16_t nucleotide) const;
		const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &HomopolymerDistribution(uint16_t nucleotide) const;

		// Main functions
		bool ReadBam( const char *bam_file, const char *adapter_file, const char *adapter_matrix, const char *variant_file, uint16_t num_threads, bool calculate_bias = true ); // Fill the class with the information from a bam file

		bool Load( const char *archive_file );
		bool Save( const char *archive_file ) const;
	};

}

#endif // DATASTATSINTERFACE_H
