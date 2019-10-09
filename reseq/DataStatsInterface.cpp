#include "DataStatsInterface.h"
using reseq::DataStatsInterface;
using reseq::DataStats;

//include <map>
using std::map;
//include <string>
using std::string;
//include <utility>
using std::pair;
//include <vector>
using std::vector;

//include "utilities.hpp"
using reseq::intQualDiff;
using reseq::uintPercent;
using reseq::uintQual;
using reseq::uintTile;
using reseq::uintSeqLen;
using reseq::uintFragCount;
using reseq::uintNucCount;

DataStatsInterface::DataStatsInterface(Reference *ref):
	stats_(ref)
	{
}

const char *DataStatsInterface::AdapterName(uintTempSeq template_segment, uintAdapterId id) const{
	return stats_.Adapters().Name(template_segment, id).c_str();
}

const vector<uintFragCount> &DataStatsInterface::AdapterCount(uintTempSeq template_segment) const{
	return stats_.Adapters().Counts(template_segment);
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::AdapterPolyATailLength() const{
	return stats_.Adapters().PolyATailLength().std();
}

uintNucCount DataStatsInterface::AdapterOverrunBases(uintBaseCall nucleotide) const{
	return stats_.Adapters().OverrunBases().at(nucleotide);
}

uintSeqLen DataStatsInterface::ErrorRatesByDistanceStart() const{
	return stats_.Coverage().ErrorRatesByDistanceSum().from();
}
uintSeqLen DataStatsInterface::ErrorRatesByDistanceEnd() const{
	return stats_.Coverage().ErrorRatesByDistanceSum().to();
}
const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &DataStatsInterface::ErrorRatesByDistance(uintSeqLen distance) const{
	return stats_.Coverage().ErrorRatesByDistanceSum()[distance].std();
}

uintPercent DataStatsInterface::ErrorRatesByGCStart() const{
	return stats_.Coverage().ErrorRatesByGCSum().from();
}
uintPercent DataStatsInterface::ErrorRatesByGCEnd() const{
	return stats_.Coverage().ErrorRatesByGCSum().to();
}
const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &DataStatsInterface::ErrorRatesByGC(uintPercent gc) const{
	return stats_.Coverage().ErrorRatesByGCSum()[gc].std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::Coverage() const{
	return stats_.Coverage().Coverage().std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::CoverageStranded( bool reverse_strand ) const{
	return stats_.Coverage().CoverageStranded(reverse_strand).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::CoverageStrandedPercent( bool reverse_strand ) const{
	return stats_.Coverage().CoverageStrandedPercent(reverse_strand).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::CoverageStrandedPercentMinCov10( bool reverse_strand ) const{
	return stats_.Coverage().CoverageStrandedPercentMinCov10(reverse_strand).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::CoverageStrandedPercentMinCov20( bool reverse_strand ) const{
	return stats_.Coverage().CoverageStrandedPercentMinCov20(reverse_strand).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::ErrorCoverage() const{
	return stats_.Coverage().ErrorCoverage().std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::ErrorCoveragePercent() const{
	return stats_.Coverage().ErrorCoveragePercent().std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::ErrorCoveragePercentMinCov10() const{
	return stats_.Coverage().ErrorCoveragePercentMinCov10().std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::ErrorCoveragePercentMinCov20() const{
	return stats_.Coverage().ErrorCoveragePercentMinCov20().std();
}

uintNucCount DataStatsInterface::ErrorCoveragePercentStranded( uintPercent forward_errors_percent, uintPercent reverse_errors_percent ) const{
	return stats_.Coverage().ErrorCoveragePercentStranded()[forward_errors_percent][reverse_errors_percent];
}

uintNucCount DataStatsInterface::ErrorCoveragePercentStrandedMinCov10( uintPercent forward_errors_percent, uintPercent reverse_errors_percent ) const{
	return stats_.Coverage().ErrorCoveragePercentStrandedMinCov10()[forward_errors_percent][reverse_errors_percent];
}

uintNucCount DataStatsInterface::ErrorCoveragePercentStrandedMinCov20( uintPercent forward_errors_percent, uintPercent reverse_errors_percent ) const{
	return stats_.Coverage().ErrorCoveragePercentStrandedMinCov20()[forward_errors_percent][reverse_errors_percent];
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::ErrorsPerRead(uintTempSeq template_segment) const{
	return stats_.Errors().ErrorsPerRead(template_segment).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::CalledBasesByBaseQualityPerPreviousCalledBase(uintTempSeq template_segment, uintBaseCall ref_base, uintBaseCall called_base, uintBaseCall previous_base) const{
	return stats_.Errors().CalledBasesByBaseQualityPerPreviousCalledBase(template_segment, ref_base, called_base, previous_base).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::CalledBasesByBaseQuality(uintTempSeq template_segment, uintBaseCall ref_base, uintBaseCall called_base) const{
	return stats_.Errors().CalledBasesByBaseQuality(template_segment, ref_base, called_base).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::CalledBasesByPosition(uintTempSeq template_segment, uintBaseCall ref_base, uintBaseCall called_base) const{
	return stats_.Errors().CalledBasesByPosition(template_segment, ref_base, called_base).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::InDelErrorByLength(uintInDelType indel_type) const{
	return stats_.Errors().InDelErrorByLength(indel_type).std();
}
const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::InDelErrorByPosition(uintInDelType indel_type) const{
	return stats_.Errors().InDelErrorByPosition(indel_type).std();
}
const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::InDelErrorByGC(uintInDelType indel_type) const{
	return stats_.Errors().InDelErrorByGC(indel_type).std();
}

const std::vector<double> &DataStatsInterface::DispersionList() const{
	return stats_.Duplicates().DispersionList();
}
const std::vector<double> &DataStatsInterface::MeanList() const{
	return stats_.Duplicates().MeanList();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::FragmentDuplicationNumber() const{
	return stats_.Duplicates().DuplicationNumber().std();
}

const vector<uintFragCount> &DataStatsInterface::Abundance() const{
	return stats_.FragmentDistribution().Abundance();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::InsertLengths() const{
	return stats_.FragmentDistribution().InsertLengths().std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::GCFragmentContent() const{
	return stats_.FragmentDistribution().GCFragmentContent().std();
}

const pair< vector<double>::size_type, vector<double> > &DataStatsInterface::GCFragmentContentBias() const{
	return stats_.FragmentDistribution().GCFragmentContentBias().std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::OutskirtContent(uintTempSeq direction, uintBaseCall nucleotide) const{
	return stats_.FragmentDistribution().OutskirtContent(direction, nucleotide).std();
}

const vector<double> &DataStatsInterface::FragmentSurroundingBiasByBase(uintBaseCall nucleotide) const{
	return stats_.FragmentDistribution().FragmentSurroundingBiasByBase(nucleotide);
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityProbabilityMean(uintTempSeq template_segment) const{
	return stats_.Qualities().SequenceQualityProbabilityMean(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityMinimum(uintTempSeq template_segment) const{
	return stats_.Qualities().SequenceQualityMinimum(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityFirstQuartile(uintTempSeq template_segment) const{
	return stats_.Qualities().SequenceQualityFirstQuartile(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityMedian(uintTempSeq template_segment) const{
	return stats_.Qualities().SequenceQualityMedian(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityThirdQuartile(uintTempSeq template_segment) const{
	return stats_.Qualities().SequenceQualityThirdQuartile(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityMaximum(uintTempSeq template_segment) const{
	return stats_.Qualities().SequenceQualityMaximum(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityContent(uintTempSeq template_segment, uintQual quality) const{
	return stats_.Qualities().SequenceQualityContent(template_segment)[quality].std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::HomoqualityDistribution(uintQual quality) const{
	return stats_.Qualities().HomoqualityDistribution()[quality].std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::NucleotideQuality(uintTempSeq template_segment, uintBaseCall nucleotide) const{
	return stats_.Qualities().NucleotideQuality(template_segment, nucleotide).stdQualities();
}

const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &DataStatsInterface::BaseQualityStatsReference(uintTempSeq template_segment, uintReadLen read_position) const{
	return stats_.Qualities().BaseQualityStatsReference(template_segment)[read_position].stdQualities();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMeanReference(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityMeanReference(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMinimumReference(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityMinimumReference(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityFirstQuartileReference(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityFirstQuartileReference(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMedianReference(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityMedianReference(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityThirdQuartileReference(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityThirdQuartileReference(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMaximumReference(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityMaximumReference(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::AverageSequenceQualityForGC(uintTempSeq template_segment) const{
	return stats_.Qualities().AverageSequenceQualityForGC(template_segment).std();
}

const std::pair< std::vector<uintNucCount>::size_type, std::vector<uintNucCount> > &DataStatsInterface::BaseQualityStats(uintTempSeq template_segment, uintReadLen read_position) const{
	return stats_.Qualities().BaseQualityStats(template_segment)[read_position].stdQualities();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMean(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityMean(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMinimum(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityMinimum(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityFirstQuartile(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityFirstQuartile(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMedian(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityMedian(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityThirdQuartile(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityThirdQuartile(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMaximum(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityMaximum(template_segment).std();
}

const pair< vector<intQualDiff>::size_type, vector<intQualDiff> > &DataStatsInterface::TileQualityMeanDifference(uintTempSeq template_segment, uintTileId tile_id) const{
	return stats_.Qualities().TileQualityMeanDifference(template_segment)[tile_id].std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::BaseQualityMeanPerStrand(uintTempSeq strand) const{
	return stats_.Qualities().BaseQualityMeanPerStrand(strand).std();
}

uintQual DataStatsInterface::BaseQualityForSequenceStart(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityForSequence(template_segment).from();
}

uintQual DataStatsInterface::BaseQualityForSequenceEnd(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityForSequence(template_segment).to();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::BaseQualityForSequence(uintTempSeq template_segment, uintQual seq_quality) const{
	return stats_.Qualities().BaseQualityForSequence(template_segment)[seq_quality].std();
}

uintQual DataStatsInterface::BaseQualityForPrecedingQualityStart(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityForPrecedingQuality(template_segment).from();
}

uintQual DataStatsInterface::BaseQualityForPrecedingQualityEnd(uintTempSeq template_segment) const{
	return stats_.Qualities().BaseQualityForPrecedingQuality(template_segment).to();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::BaseQualityForPrecedingQuality(uintTempSeq template_segment, uintQual seq_quality) const{
	return stats_.Qualities().BaseQualityForPrecedingQuality(template_segment)[seq_quality].std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityMean(uintTempSeq template_segment) const{
	return stats_.Qualities().SequenceQualityMean(template_segment).std();
}

uintQual DataStatsInterface::SequenceQualityPairsStart() const{
	return stats_.Qualities().SequenceQualityMeanPaired().from();
}

uintQual DataStatsInterface::SequenceQualityPairsEnd() const{
	return stats_.Qualities().SequenceQualityMeanPaired().to();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SequenceQualityPairs(uintQual seq_quality_first) const{
	return stats_.Qualities().SequenceQualityMeanPaired()[seq_quality_first].std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::MeanSequenceQualityMeanByFragmentLength(uintTempSeq template_segment) const{
	return stats_.Qualities().MeanSequenceQualityMeanByFragmentLength(template_segment).std();
}

const pair< vector<uintQual>::size_type, vector<uintQual> > &DataStatsInterface::AverageSequenceQualityForBase(uintTempSeq template_segment, uintBaseCall nucleotide) const{
	return stats_.Qualities().AverageSequenceQualityForBase(template_segment, nucleotide).std();
}

const vector<uintTile> &DataStatsInterface::TileNames() const{
	return stats_.Tiles().Tiles();
}

const vector<uintFragCount> &DataStatsInterface::TileAbundance() const{
	return stats_.Tiles().Abundance();
}

uintQual DataStatsInterface::PhredQualityOffset() const{
	return stats_.PhredQualityOffset();
}
uintFragCount DataStatsInterface::TotalNumberReads() const{
	return stats_.TotalNumberReads();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::ReadLengths(uintTempSeq template_segment) const{
	return stats_.ReadLengths(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::ProperPairMappingQuality() const{
	return stats_.ProperPairMappingQuality().std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::ImproperPairMappingQuality() const{
	return stats_.ImproperPairMappingQuality().std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::SingleReadMappingQuality() const{
	return stats_.SingleReadMappingQuality().std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::GCReadContent(uintTempSeq template_segment) const{
	return stats_.GCReadContent(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::GCReadContentReference(uintTempSeq template_segment) const{
	return stats_.GCReadContentReference(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::GCReadContentMapped(uintTempSeq template_segment) const{
	return stats_.GCReadContentMapped(template_segment).std();
}

const pair< vector<uintFragCount>::size_type, vector<uintFragCount> > &DataStatsInterface::NContent(uintTempSeq template_segment) const{
	return stats_.NContent(template_segment).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::SequenceContent(uintTempSeq template_segment, uintBaseCall nucleotide) const{
	return stats_.SequenceContent(template_segment, nucleotide).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::SequenceContentReference(uintTempSeq template_segment, uintTempSeq strand, uintBaseCall nucleotide) const{
	return stats_.SequenceContentReference(template_segment, strand, nucleotide).std();
}

const pair< vector<uintNucCount>::size_type, vector<uintNucCount> > &DataStatsInterface::HomopolymerDistribution(uintBaseCall nucleotide) const{
	return stats_.HomopolymerDistribution(nucleotide).std();
}

// Deactivation of pcr calculation only for speeding up of tests
bool DataStatsInterface::ReadBam( const char *bam_file, const char *adapter_file, const char *adapter_matrix, const char *variant_file, uintNumThreads num_threads, bool calculate_bias ){
	bool success = stats_.ReadBam( bam_file, adapter_file, adapter_matrix, string(variant_file), num_threads, calculate_bias );
	stats_.PreparePlotting();
	return success;
}

bool DataStatsInterface::Load( const char *archive_file ){
	bool success = stats_.Load( archive_file );
	stats_.PreparePlotting();
	return success;
}

bool DataStatsInterface::Save( const char *archive_file ) const{
	return stats_.Save( archive_file );
}
