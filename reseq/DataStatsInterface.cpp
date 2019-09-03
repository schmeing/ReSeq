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

DataStatsInterface::DataStatsInterface(Reference *ref):
	stats_(ref)
	{
}

const char *DataStatsInterface::AdapterName(uint16_t template_segment, uint16_t id) const{
	return stats_.Adapters().Name(template_segment, id).c_str();
}

const vector<uint64_t> &DataStatsInterface::AdapterCount(uint16_t template_segment) const{
	return stats_.Adapters().Counts(template_segment);
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::AdapterPolyATailLength() const{
	return stats_.Adapters().PolyATailLength().std();
}

uint64_t DataStatsInterface::AdapterOverrunBases(uint16_t nucleotide) const{
	return stats_.Adapters().OverrunBases().at(nucleotide);
}

uint32_t DataStatsInterface::ErrorRatesByDistanceStart() const{
	return stats_.Coverage().ErrorRatesByDistanceSum().from();
}
uint32_t DataStatsInterface::ErrorRatesByDistanceEnd() const{
	return stats_.Coverage().ErrorRatesByDistanceSum().to();
}
const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &DataStatsInterface::ErrorRatesByDistance(uint32_t distance) const{
	return stats_.Coverage().ErrorRatesByDistanceSum()[distance].std();
}

uint16_t DataStatsInterface::ErrorRatesByGCStart() const{
	return stats_.Coverage().ErrorRatesByGCSum().from();
}
uint16_t DataStatsInterface::ErrorRatesByGCEnd() const{
	return stats_.Coverage().ErrorRatesByGCSum().to();
}
const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &DataStatsInterface::ErrorRatesByGC(uint16_t gc) const{
	return stats_.Coverage().ErrorRatesByGCSum()[gc].std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::Coverage() const{
	return stats_.Coverage().Coverage().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::CoverageStranded( bool reverse_strand ) const{
	return stats_.Coverage().CoverageStranded(reverse_strand).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::CoverageStrandedPercent( bool reverse_strand ) const{
	return stats_.Coverage().CoverageStrandedPercent(reverse_strand).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::CoverageStrandedPercentMinCov10( bool reverse_strand ) const{
	return stats_.Coverage().CoverageStrandedPercentMinCov10(reverse_strand).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::CoverageStrandedPercentMinCov20( bool reverse_strand ) const{
	return stats_.Coverage().CoverageStrandedPercentMinCov20(reverse_strand).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::ErrorCoverage() const{
	return stats_.Coverage().ErrorCoverage().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::ErrorCoveragePercent() const{
	return stats_.Coverage().ErrorCoveragePercent().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::ErrorCoveragePercentMinCov10() const{
	return stats_.Coverage().ErrorCoveragePercentMinCov10().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::ErrorCoveragePercentMinCov20() const{
	return stats_.Coverage().ErrorCoveragePercentMinCov20().std();
}

uint64_t DataStatsInterface::ErrorCoveragePercentStranded( unsigned char forward_errors_percent, unsigned char reverse_errors_percent ) const{
	return stats_.Coverage().ErrorCoveragePercentStranded()[forward_errors_percent][reverse_errors_percent];
}

uint64_t DataStatsInterface::ErrorCoveragePercentStrandedMinCov10( unsigned char forward_errors_percent, unsigned char reverse_errors_percent ) const{
	return stats_.Coverage().ErrorCoveragePercentStrandedMinCov10()[forward_errors_percent][reverse_errors_percent];
}

uint64_t DataStatsInterface::ErrorCoveragePercentStrandedMinCov20( unsigned char forward_errors_percent, unsigned char reverse_errors_percent ) const{
	return stats_.Coverage().ErrorCoveragePercentStrandedMinCov20()[forward_errors_percent][reverse_errors_percent];
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::ErrorsPerRead(uint16_t template_segment) const{
	return stats_.Errors().ErrorsPerRead(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::CalledBasesByBaseQualityPerPreviousCalledBase(uint16_t template_segment, uint16_t ref_base, uint16_t called_base, uint16_t previous_base) const{
	return stats_.Errors().CalledBasesByBaseQualityPerPreviousCalledBase(template_segment, ref_base, called_base, previous_base).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::CalledBasesByBaseQuality(uint16_t template_segment, uint16_t ref_base, uint16_t called_base) const{
	return stats_.Errors().CalledBasesByBaseQuality(template_segment, ref_base, called_base).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::CalledBasesByPosition(uint16_t template_segment, uint16_t ref_base, uint16_t called_base) const{
	return stats_.Errors().CalledBasesByPosition(template_segment, ref_base, called_base).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::InDelErrorByLength(uint16_t indel_type) const{
	return stats_.Errors().InDelErrorByLength(indel_type).std();
}
const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::InDelErrorByPosition(uint16_t indel_type) const{
	return stats_.Errors().InDelErrorByPosition(indel_type).std();
}
const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::InDelErrorByGC(uint16_t indel_type) const{
	return stats_.Errors().InDelErrorByGC(indel_type).std();
}

const std::vector<double> &DataStatsInterface::DispersionList() const{
	return stats_.Duplicates().DispersionList();
}
const std::vector<double> &DataStatsInterface::MeanList() const{
	return stats_.Duplicates().MeanList();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::FragmentDuplicationNumber() const{
	return stats_.Duplicates().DuplicationNumber().std();
}

const vector<uint64_t> &DataStatsInterface::Abundance() const{
	return stats_.FragmentDistribution().Abundance();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::InsertLengths() const{
	return stats_.FragmentDistribution().InsertLengths().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::GCFragmentContent() const{
	return stats_.FragmentDistribution().GCFragmentContent().std();
}

const pair< vector<double>::size_type, vector<double> > &DataStatsInterface::GCFragmentContentBias() const{
	return stats_.FragmentDistribution().GCFragmentContentBias().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::OutskirtContent(uint16_t direction, uint16_t nucleotide) const{
	return stats_.FragmentDistribution().OutskirtContent(direction, nucleotide).std();
}

const vector<double> &DataStatsInterface::FragmentSurroundingBiasByBase(uint16_t nucleotide) const{
	return stats_.FragmentDistribution().FragmentSurroundingBiasByBase(nucleotide);
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityProbabilityMean(uint16_t template_segment) const{
	return stats_.Qualities().SequenceQualityProbabilityMean(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityMinimum(uint16_t template_segment) const{
	return stats_.Qualities().SequenceQualityMinimum(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityFirstQuartile(uint16_t template_segment) const{
	return stats_.Qualities().SequenceQualityFirstQuartile(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityMedian(uint16_t template_segment) const{
	return stats_.Qualities().SequenceQualityMedian(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityThirdQuartile(uint16_t template_segment) const{
	return stats_.Qualities().SequenceQualityThirdQuartile(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityMaximum(uint16_t template_segment) const{
	return stats_.Qualities().SequenceQualityMaximum(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityContent(uint16_t template_segment, unsigned char quality) const{
	return stats_.Qualities().SequenceQualityContent(template_segment)[quality].std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::HomoqualityDistribution(uint16_t quality) const{
	return stats_.Qualities().HomoqualityDistribution()[quality].std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::NucleotideQuality(uint16_t template_segment, uint16_t nucleotide) const{
	return stats_.Qualities().NucleotideQuality(template_segment, nucleotide).stdQualities();
}

const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &DataStatsInterface::BaseQualityStatsReference(uint16_t template_segment, uint32_t read_position) const{
	return stats_.Qualities().BaseQualityStatsReference(template_segment)[read_position].stdQualities();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMeanReference(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityMeanReference(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMinimumReference(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityMinimumReference(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityFirstQuartileReference(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityFirstQuartileReference(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMedianReference(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityMedianReference(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityThirdQuartileReference(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityThirdQuartileReference(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMaximumReference(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityMaximumReference(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::AverageSequenceQualityForGC(uint16_t template_segment) const{
	return stats_.Qualities().AverageSequenceQualityForGC(template_segment).std();
}

const std::pair< std::vector<uint64_t>::size_type, std::vector<uint64_t> > &DataStatsInterface::BaseQualityStats(uint16_t template_segment, uint32_t read_position) const{
	return stats_.Qualities().BaseQualityStats(template_segment)[read_position].stdQualities();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMean(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityMean(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMinimum(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityMinimum(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityFirstQuartile(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityFirstQuartile(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMedian(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityMedian(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityThirdQuartile(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityThirdQuartile(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMaximum(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityMaximum(template_segment).std();
}

const pair< vector<signed char>::size_type, vector<signed char> > &DataStatsInterface::TileQualityMeanDifference(uint16_t template_segment, uint16_t tile_id) const{
	return stats_.Qualities().TileQualityMeanDifference(template_segment)[tile_id].std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::BaseQualityMeanPerStrand(uint16_t strand) const{
	return stats_.Qualities().BaseQualityMeanPerStrand(strand).std();
}

uint16_t DataStatsInterface::BaseQualityForSequenceStart(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityForSequence(template_segment).from();
}

uint16_t DataStatsInterface::BaseQualityForSequenceEnd(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityForSequence(template_segment).to();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::BaseQualityForSequence(uint16_t template_segment, uint16_t seq_quality) const{
	return stats_.Qualities().BaseQualityForSequence(template_segment)[seq_quality].std();
}

uint16_t DataStatsInterface::BaseQualityForPrecedingQualityStart(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityForPrecedingQuality(template_segment).from();
}

uint16_t DataStatsInterface::BaseQualityForPrecedingQualityEnd(uint16_t template_segment) const{
	return stats_.Qualities().BaseQualityForPrecedingQuality(template_segment).to();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::BaseQualityForPrecedingQuality(uint16_t template_segment, uint16_t seq_quality) const{
	return stats_.Qualities().BaseQualityForPrecedingQuality(template_segment)[seq_quality].std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityMean(uint16_t template_segment) const{
	return stats_.Qualities().SequenceQualityMean(template_segment).std();
}

uint16_t DataStatsInterface::SequenceQualityPairsStart() const{
	return stats_.Qualities().SequenceQualityMeanPaired().from();
}

uint16_t DataStatsInterface::SequenceQualityPairsEnd() const{
	return stats_.Qualities().SequenceQualityMeanPaired().to();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceQualityPairs(uint16_t seq_quality_first) const{
	return stats_.Qualities().SequenceQualityMeanPaired()[seq_quality_first].std();
}

const pair< vector<uint16_t>::size_type, vector<uint16_t> > &DataStatsInterface::MeanSequenceQualityMeanByFragmentLength(uint16_t template_segment) const{
	return stats_.Qualities().MeanSequenceQualityMeanByFragmentLength(template_segment).std();
}

const pair< vector<unsigned char>::size_type, vector<unsigned char> > &DataStatsInterface::AverageSequenceQualityForBase(uint16_t template_segment, uint16_t nucleotide) const{
	return stats_.Qualities().AverageSequenceQualityForBase(template_segment, nucleotide).std();
}

const vector<uint16_t> &DataStatsInterface::TileNames() const{
	return stats_.Tiles().Tiles();
}

const vector<uint64_t> &DataStatsInterface::TileAbundance() const{
	return stats_.Tiles().Abundance();
}

unsigned char DataStatsInterface::PhredQualityOffset() const{
	return stats_.PhredQualityOffset();
}
uint64_t DataStatsInterface::TotalNumberReads() const{
	return stats_.TotalNumberReads();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::ReadLengths(uint16_t template_segment) const{
	return stats_.ReadLengths(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::ProperPairMappingQuality() const{
	return stats_.ProperPairMappingQuality().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::ImproperPairMappingQuality() const{
	return stats_.ImproperPairMappingQuality().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SingleReadMappingQuality() const{
	return stats_.SingleReadMappingQuality().std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::GCReadContent(uint16_t template_segment) const{
	return stats_.GCReadContent(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::GCReadContentReference(uint16_t template_segment) const{
	return stats_.GCReadContentReference(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::GCReadContentMapped(uint16_t template_segment) const{
	return stats_.GCReadContentMapped(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::NContent(uint16_t template_segment) const{
	return stats_.NContent(template_segment).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceContent(uint16_t template_segment, uint16_t nucleotide) const{
	return stats_.SequenceContent(template_segment, nucleotide).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::SequenceContentReference(uint16_t template_segment, uint16_t strand, uint16_t nucleotide) const{
	return stats_.SequenceContentReference(template_segment, strand, nucleotide).std();
}

const pair< vector<uint64_t>::size_type, vector<uint64_t> > &DataStatsInterface::HomopolymerDistribution(uint16_t nucleotide) const{
	return stats_.HomopolymerDistribution(nucleotide).std();
}

// Deactivation of pcr calculation only for speeding up of tests
bool DataStatsInterface::ReadBam( const char *bam_file, const char *adapter_file, const char *adapter_matrix, const char *variant_file, uint16_t num_threads, bool calculate_bias ){
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
