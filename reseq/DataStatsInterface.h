#ifndef DATASTATSINTERFACE_H
#define DATASTATSINTERFACE_H

#include <map>
#include <stdint.h>
#include <string>
#include <utility>
#include <vector>

#include "DataStats.h"
#include "utilities.hpp"

namespace reseq{
	class DataStatsInterface{
	private:
		DataStats stats_;

	public:
		DataStatsInterface(Reference *ref);

		// Getter functions with std return values (used for python plotting)
		const char *AdapterName(uintTempSeq template_segment, uintAdapterId id) const;
		const std::vector<uintFragCount> &AdapterCount(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &AdapterPolyATailLength() const;
		uintNucCount AdapterOverrunBases(uintBaseCall nucleotide) const;

		uintSeqLen ErrorRatesByDistanceStart() const;
		uintSeqLen ErrorRatesByDistanceEnd() const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &ErrorRatesByDistance(uintSeqLen distance) const;
		uintPercent ErrorRatesByGCStart() const;
		uintPercent ErrorRatesByGCEnd() const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &ErrorRatesByGC(uintPercent gc) const;

		const std::pair< std::vector<uintSurBlockId>::size_type, std::vector<uintSurBlockId> > &BlockErrorRate() const;
		const std::pair< std::vector<uintSurBlockId>::size_type, std::vector<uintSurBlockId> > &BlockPercentSystematic() const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &SystematicErrorPValues() const;

		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &Coverage() const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &CoverageStranded( bool reverse_strand ) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &CoverageStrandedPercent( bool reverse_strand ) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &CoverageStrandedPercentMinCov10( bool reverse_strand ) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &CoverageStrandedPercentMinCov20( bool reverse_strand ) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &ErrorCoverage() const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &ErrorCoveragePercent() const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &ErrorCoveragePercentMinCov10() const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &ErrorCoveragePercentMinCov20() const;
		uintNucCount ErrorCoveragePercentStranded( uintPercent forward_errors_percent, uintPercent reverse_errors_percent ) const;
		uintNucCount ErrorCoveragePercentStrandedMinCov10( uintPercent forward_errors_percent, uintPercent reverse_errors_percent ) const;
		uintNucCount ErrorCoveragePercentStrandedMinCov20( uintPercent forward_errors_percent, uintPercent reverse_errors_percent ) const;

		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &ErrorsPerRead(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &CalledBasesByBaseQualityPerPreviousCalledBase(uintTempSeq template_segment, uintBaseCall ref_base, uintBaseCall called_base, uintBaseCall previous_base) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &CalledBasesByBaseQuality(uintTempSeq template_segment, uintBaseCall ref_base, uintBaseCall called_base) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &CalledBasesByPosition(uintTempSeq template_segment,uintBaseCall ref_base, uintBaseCall called_base) const;

		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &InDelErrorByLength(uintInDelType indel_type) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &InDelErrorByPosition(uintInDelType indel_type) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &InDelErrorByGC(uintInDelType indel_type) const;

		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &FragmentDuplicationNumber() const;

		const std::vector<uintFragCount> &Abundance() const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &InsertLengths() const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &GCFragmentContent() const;
		const std::pair< std::vector<double>::size_type, std::vector<double> > &GCFragmentContentBias() const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &OutskirtContent(uintTempSeq direction, uintBaseCall nucleotide) const;
		const std::vector<double> &FragmentSurroundingBiasByBase(uintBaseCall nucleotide) const;

		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityProbabilityMean(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityMinimum(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityFirstQuartile(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityMedian(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityThirdQuartile(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityMaximum(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityContent(uintTempSeq template_segment, uintQual quality) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &HomoqualityDistribution(uintQual quality) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &NucleotideQuality(uintTempSeq template_segment, uintBaseCall nucleotide) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &BaseQualityStatsReference(uintTempSeq template_segment, uintReadLen read_position) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMeanReference(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMinimumReference(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityFirstQuartileReference(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMedianReference(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityThirdQuartileReference(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMaximumReference(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &AverageSequenceQualityForGC(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &BaseQualityStats(uintTempSeq template_segment, uintReadLen read_position) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMean(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMinimum(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityFirstQuartile(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMedian(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityThirdQuartile(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMaximum(uintTempSeq template_segment) const;
		const std::pair< std::vector<intQualDiff>::size_type, std::vector<intQualDiff> > &TileQualityMeanDifference(uintTempSeq template_segment, uintTileId tile_id) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &BaseQualityMeanPerStrand(uintTempSeq strand) const;
		uintQual BaseQualityForSequenceStart(uintTempSeq template_segment) const;
		uintQual BaseQualityForSequenceEnd(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &BaseQualityForSequence(uintTempSeq template_segment, uintQual seq_quality) const;
		uintQual BaseQualityForPrecedingQualityStart(uintTempSeq template_segment) const;
		uintQual BaseQualityForPrecedingQualityEnd(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &BaseQualityForPrecedingQuality(uintTempSeq template_segment, uintQual seq_quality) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityMean(uintTempSeq template_segment) const;
		uintQual SequenceQualityPairsStart() const;
		uintQual SequenceQualityPairsEnd() const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SequenceQualityPairs(uintQual seq_quality_first) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &MeanSequenceQualityMeanByFragmentLength(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintQual>::size_type, std::vector<uintQual> > &AverageSequenceQualityForBase(uintTempSeq template_segment, uintBaseCall nucleotide) const;

		const std::vector<uintTile> &TileNames() const;
		const std::vector<uintFragCount> &TileAbundance() const;

		uintQual PhredQualityOffset() const;
		uintFragCount TotalNumberReads() const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &ReadLengths(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &ProperPairMappingQuality() const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &ImproperPairMappingQuality() const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &SingleReadMappingQuality() const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &GCReadContent(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &GCReadContentReference(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &GCReadContentMapped(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintFragCount>::size_type, std::vector<uintFragCount> > &NContent(uintTempSeq template_segment) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &SequenceContent(uintTempSeq template_segment, uintBaseCall nucleotide) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &SequenceContentReference(uintTempSeq template_segment, uintTempSeq strand, uintBaseCall nucleotide) const;
		const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &HomopolymerDistribution(uintBaseCall nucleotide) const;

		// Main functions
		bool ReadBam( const char *bam_file, const char *adapter_file, const char *adapter_matrix, const char *variant_file, uintNumThreads num_threads, bool calculate_bias = true ); // Fill the class with the information from a bam file

		bool Load( const char *archive_file );
		bool Save( const char *archive_file ) const;
	};

}

#endif // DATASTATSINTERFACE_H
