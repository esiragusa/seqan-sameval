// ==========================================================================
//                                SAMeval
// ==========================================================================
// Copyright (C) 2014 Enrico Siragusa, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

struct Options
{
    CharString inputOracle;
    CharString inputMapper;
    CharString outputTsv;

    bool precision;
    bool verbose;
    bool delta;

    Options() :
        precision(false),
        verbose(false),
        delta(false)
    {}
};

// ----------------------------------------------------------------------------
// Class RabemaOracleStats
// ----------------------------------------------------------------------------

namespace seqan {
struct Tsv_;
typedef Tag<Tsv_> Tsv;
}

// ----------------------------------------------------------------------------
// Class RabemaStats
// ----------------------------------------------------------------------------

struct RabemaStats
{
    // Number of intervals that were to find.
    __uint64 intervalsToFind;
    // Number of intervals that were actually found.
    __uint64 intervalsFound;
    // SAM records that did not correspond to alignments below the configured maximal error rate.
    __uint64 invalidAlignments;

    // Total number of reads.
    __uint64 totalReads;
    // Number of reads with alignments below maximal error rate.
    __uint64 mappedReads;
    // Total number of reads with GSI records equals number of normalized intervals to find.
    __uint64 readsInGsi;
    // Normalized number of found intervals.
    double normalizedIntervals;
    // Number of additional alignments in SAM file with low enough error rate but no GSI record.
    unsigned additionalHits;

    // The following arrays are indexed by the integer value of the error rate.

    // Number of intervals that were to find for each error rate.
    String<unsigned> intervalsToFindForErrorRate;
    // Number of found intervals for each error rate.
    String<unsigned> intervalsFoundForErrorRate;

    // The following values are normalized towards all intervals.

    // Normalized number of intervals to find for each error rate.
    String<double> normalizedIntervalsToFindForErrorRate;
    // Normalized number of intervals found for each error rate.
    String<double> normalizedIntervalsFoundForErrorRate;

    RabemaStats() :
        intervalsToFind(0), intervalsFound(0), invalidAlignments(0), totalReads(0), mappedReads(0),
        readsInGsi(0), normalizedIntervals(0), additionalHits(0)
    {}

    RabemaStats(unsigned maxErrorRate) :
        intervalsToFind(0), intervalsFound(0), invalidAlignments(0), totalReads(0), mappedReads(0),
        readsInGsi(0), normalizedIntervals(0),
        additionalHits(0)
    {
        resize(intervalsToFindForErrorRate, maxErrorRate + 1, 0);
        resize(intervalsFoundForErrorRate, maxErrorRate + 1, 0);
        resize(normalizedIntervalsToFindForErrorRate, maxErrorRate + 1, 0.0);
        resize(normalizedIntervalsFoundForErrorRate, maxErrorRate + 1, 0.0);
    }

};

struct RabemaOracleStats : RabemaStats
{
    String<unsigned> readsMappedPerErrorRate;

    RabemaOracleStats() : RabemaStats() {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser, Options const & /* options */)
{
    setAppName(parser, "rabema_oracle");
    setShortDescription(parser, "RABEMA Oracle");
    setCategory(parser, "Benchmarking");

//    setDateAndVersion(parser);
//    setDescription(parser);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIORACLE.{sam|bam}\\fP \\fIMAPPING.{sam|bam}\\fP");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "ORACLE"));
    setValidValues(parser, 0, BamFileIn::getFileExtensions());
    setHelpText(parser, 0, "Oracle SAM/BAM file to use as gold standard.");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "MAPPING"));
    setValidValues(parser, 1, BamFileIn::getFileExtensions());
    setHelpText(parser, 1, "Read mapper SAM/BAM file to evaluate.");

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("p", "precision", "Compute precision. Default: compute recall."));
    addOption(parser, seqan::ArgParseOption("d", "delta", "Tolerate +- delta bp for begin/end positions. \
                                                           Default: require exact begin/end positions."));

    addSection(parser, "Input / Output");

    addOption(parser, seqan::ArgParseOption("", "out-tsv", "Write the statistics to this TSV file.",
                                            seqan::ArgParseArgument::OUTPUT_FILE, "TSV"));
    setValidValues(parser, "out-tsv", "rabema_report_tsv");
}

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

ArgumentParser::ParseResult
parseCommandLine(Options & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.inputOracle, parser, 0);
    getArgumentValue(options.inputMapper, parser, 1);

    getOptionValue(options.verbose, parser, "verbose");
    getOptionValue(options.precision, parser, "precision");
    getOptionValue(options.delta, parser, "delta");
    getOptionValue(options.outputTsv, parser, "out-tsv");

    return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function getEditDistance()
// ----------------------------------------------------------------------------

unsigned getEditDistance(BamAlignmentRecord const & record)
{
    BamTagsDict bamTags(const_cast<CharString &>(record.tags));

    unsigned idx = 0;
    unsigned editDistance = 0;

    if (findTagKey(idx, bamTags, "NM") && extractTagValue(editDistance, bamTags, idx))
        return editDistance;
    else
        throw RuntimeError("Oracle SAM/BAM record does not contain the NM tag.");
}

// ----------------------------------------------------------------------------
// Function getErrorRate()
// ----------------------------------------------------------------------------

unsigned getErrorRate(BamAlignmentRecord const & record)
{
    return getEditDistance(record) * 100.0 / length(record.seq);
}

// ----------------------------------------------------------------------------
// Function getEndPos()
// ----------------------------------------------------------------------------

unsigned getEndPos(BamAlignmentRecord const & record)
{
    unsigned mapperRecordLength = getAlignmentLengthInRef(record) - countPaddings(record.cigar);
    return hasFlagRC(record) ? record.beginPos + mapperRecordLength : record.beginPos - mapperRecordLength;
}

// ----------------------------------------------------------------------------
// Function getRange()
// ----------------------------------------------------------------------------

template <typename TValue, typename TDelta>
Pair<TValue> getRange(TValue value, TDelta delta)
{
    SEQAN_ASSERT_GEQ(value, delta);
    return Pair<TValue>(value - delta, value + delta);
}

// ----------------------------------------------------------------------------
// Function isInRange()
// ----------------------------------------------------------------------------

template <typename TValue, typename TRange>
bool isInRange(TValue value, Pair<TRange> range)
{
    return value >= range.i1 && value <= range.i2;
}

// ----------------------------------------------------------------------------
// Function isEqual()
// ----------------------------------------------------------------------------
// Check contig id, orientation and position.

bool isEqual(BamAlignmentRecord const & mapperRecord, BamAlignmentRecord const & oracleRecord, unsigned delta = 0)
{
    return (mapperRecord.rID == oracleRecord.rID) &&
           (hasFlagRC(mapperRecord) == hasFlagRC(oracleRecord)) &&
           (isInRange(mapperRecord.beginPos, getRange(oracleRecord.beginPos, delta))) &&
           (isInRange(getEndPos(mapperRecord), getRange(getEndPos(oracleRecord), delta)));
}

// ----------------------------------------------------------------------------
// Function updateMaximalErrorRate()
// ----------------------------------------------------------------------------

void updateMaximalErrorRate(RabemaStats & stats, unsigned maxErrorRate)
{
    if (length(stats.intervalsToFindForErrorRate) <= maxErrorRate + 1)
        resize(stats.intervalsToFindForErrorRate, maxErrorRate + 1, 0);
    if (length(stats.intervalsFoundForErrorRate) <= maxErrorRate + 1)
        resize(stats.intervalsFoundForErrorRate, maxErrorRate + 1, 0);
    if (length(stats.normalizedIntervalsToFindForErrorRate) <= maxErrorRate + 1)
        resize(stats.normalizedIntervalsToFindForErrorRate, maxErrorRate + 1, 0.0);
    if (length(stats.normalizedIntervalsFoundForErrorRate) <= maxErrorRate + 1)
        resize(stats.normalizedIntervalsFoundForErrorRate, maxErrorRate + 1, 0.0);
}

void updateMaximalErrorRate(RabemaOracleStats & stats, unsigned maxErrorRate)
{
    updateMaximalErrorRate(static_cast<RabemaStats &>(stats), maxErrorRate);

    if (length(stats.readsMappedPerErrorRate) <= maxErrorRate + 1)
        resize(stats.readsMappedPerErrorRate, maxErrorRate + 1, 0);
}

// ----------------------------------------------------------------------------
// Function writeStats()
// ----------------------------------------------------------------------------

template <typename TStream>
void writeStats(TStream & stream, RabemaStats const & stats)
{
    stream << "Intervals to find:              " << stats.intervalsToFind << '\n'
           << "Intervals found:                " << stats.intervalsFound << '\n';
    if (stats.intervalsToFind > 0u)
        stream << "Intervals found [%]             " << (100.0 * stats.intervalsFound / stats.intervalsToFind) << '\n';
    else
        stream << "Intervals found [%]             " << 0 << '\n';
    stream << "Invalid alignments:             " << stats.invalidAlignments << '\n'
           << "Additional Hits:                " << stats.additionalHits << '\n'
           << '\n'
           << "Number of reads:                " << stats.totalReads << '\n'
           << "Number of reads with intervals: " << stats.readsInGsi << '\n'
           << '\n'
           << "Mapped reads:                   " << stats.mappedReads << '\n'
           << "Mapped reads [% of total]:      " << 100.0 * stats.mappedReads / stats.totalReads << '\n'
           << "Mapped reads [% of mappable]:   " << 100.0 * stats.mappedReads / stats.readsInGsi  << '\n'
           << '\n'
           << "Normalized intervals found:     " << stats.normalizedIntervals << '\n';
    if (stats.readsInGsi > 0u)
        stream << "Normalized intervals found [%]: " << (100.0 * stats.normalizedIntervals / stats.readsInGsi) << "\n\n";
    else
        stream << "Normalized intervals found [%]: " << 0 << "\n\n";
    stream << "Found Intervals By Error Rate\n"
           << "\n";
    char buffer[1000];
    sprintf(buffer, "  ERR\t%8s\t%8s\t%8s\t%8s\t%10s\t%10s\n", "#max", "#found", "%found", "norm max", "norm found", "norm found [%]");
    stream << buffer
           << "------------------------------------------------------------------------------------------------------\n";
    for (unsigned i = 0; i < length(stats.intervalsToFindForErrorRate); ++i)
    {
        // if (maxError != -1  && (int)i != maxError)
        //     continue;
        double percFoundIntervals = 100.0 * stats.intervalsFoundForErrorRate[i] / stats.intervalsToFindForErrorRate[i];
        if (stats.intervalsToFindForErrorRate[i] == 0u)
            percFoundIntervals = 0;
        double percFoundNormalizedIntervals = 100.0 * stats.normalizedIntervalsFoundForErrorRate[i] / stats.normalizedIntervalsToFindForErrorRate[i];
        if (stats.normalizedIntervalsToFindForErrorRate[i] == 0)
            percFoundNormalizedIntervals = 0;
        sprintf(buffer, "%5u\t%8d\t%8d\t%8.2f\t%8.2f\t%10.2f\t%10.2f\n", i, stats.intervalsToFindForErrorRate[i], stats.intervalsFoundForErrorRate[i],
                percFoundIntervals, stats.normalizedIntervalsToFindForErrorRate[i], stats.normalizedIntervalsFoundForErrorRate[i],
                percFoundNormalizedIntervals);
        stream << buffer;
    }
    stream << '\n';
}

// ----------------------------------------------------------------------------
// Function writeTsv()
// ----------------------------------------------------------------------------

template <typename TStream, typename TString1, typename TString2>
int writeTsv(TStream & stream, RabemaStats const & stats, int maxError, TString1 const & benchmarkCategory,
              bool oracleMode, TString2 const & distanceMetric)
{
    stream << "##Rabema Results\n"
           << "##\n"
           << "##category\t" << benchmarkCategory << '\n'
           << "##max_distance\t" << maxError<< '\n'
           << "##oracle_mode\t" << (oracleMode ? (char const *) "yes" : (char const *) "no") << "\n"
           << "##distance_metric\t" << distanceMetric << '\n'
           << "##\n"
           << "##intervals_to_find\t" << stats.intervalsToFind << '\n'
           << "##intervals_found\t" << stats.intervalsFound << '\n'
           << "##intervals_found_percent\t" << (100.0 * stats.intervalsFound / stats.intervalsToFind) << '\n'
           << "##invalid_alignments\t" << stats.invalidAlignments << '\n'
           << "##additional_hits\t" << stats.additionalHits << '\n'
           << "##\n"
           << "##number_of_reads\t" << stats.totalReads << '\n'
           << "##number_of_reads_with_intervals\t" << stats.readsInGsi << '\n'
           << "##\n"
           << "##mapped_reads\t" << stats.mappedReads << '\n'
           << "##mapped_reads_percent_of_total\t" << 100.0 * stats.mappedReads / stats.totalReads << '\n'
           << "##mapped_reads_percent_of_mappable\t" << 100.0 * stats.mappedReads / stats.readsInGsi << '\n'
           << "##\n"
           << "##normalized_intervals_found\t" << stats.normalizedIntervals << '\n'
           << "##normalized_intervals_found_percent\t" << (100.0 * stats.normalizedIntervals / stats.readsInGsi) << '\n'
           << "##\n"
           << "##Found Intervals By Distance\n"
           << "##\n"
           << "#error_rate\tnum_max\tnum_found\tpercent_found\tnorm_max\tnorm_found\tpercent_norm_found\n";
    char buffer[1000];
    for (unsigned i = 0; i < length(stats.intervalsToFindForErrorRate); ++i)
    {
        // if (maxError != -1  && (int)i != maxError)
        //     continue;
        double percFoundIntervals = 100.0 * stats.intervalsFoundForErrorRate[i] / stats.intervalsToFindForErrorRate[i];
        if (stats.intervalsToFindForErrorRate[i] == 0u)
            percFoundIntervals = 0;
        double percFoundNormalizedIntervals = 100.0 * stats.normalizedIntervalsFoundForErrorRate[i] / stats.normalizedIntervalsToFindForErrorRate[i];
        if (stats.normalizedIntervalsToFindForErrorRate[i] == 0)
            percFoundNormalizedIntervals = 0;
        sprintf(buffer, "%u\t%d\t%d\t%.2f\t%.2f\t%1.2f\t%.2f\n", i, stats.intervalsToFindForErrorRate[i], stats.intervalsFoundForErrorRate[i],
                percFoundIntervals, stats.normalizedIntervalsToFindForErrorRate[i], stats.normalizedIntervalsFoundForErrorRate[i],
                percFoundNormalizedIntervals);
        stream << buffer;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function lessThanSamtoolsQueryName()
// ----------------------------------------------------------------------------
// Original comparison function for strings by Heng Li from samtools.

inline int strnum_cmp(const char *a, const char *b)
{
	char *pa, *pb;
	pa = (char*)a; pb = (char*)b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			long ai, bi;
			ai = strtol(pa, &pa, 10);
			bi = strtol(pb, &pb, 10);
			if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
		} else {
			if (*pa != *pb) break;
			++pa; ++pb;
		}
	}
	if (*pa == *pb)
		return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
	return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}

inline bool lessThanSamtoolsQueryName(CharString const & lhs, CharString const & rhs)
{
    return strnum_cmp(toCString(lhs), toCString(rhs)) < 0;
}

// ----------------------------------------------------------------------------
// Function updateOracleStats()
// ----------------------------------------------------------------------------
// Update stats after reading one oracle record.

void updateOracleStats(RabemaOracleStats & stats, BamAlignmentRecord const & record)
{
    unsigned errorRate = getErrorRate(record);

    updateMaximalErrorRate(stats, errorRate);

    stats.totalReads++;
    stats.readsInGsi++;
    stats.intervalsToFind++;
    stats.intervalsToFindForErrorRate[errorRate]++;
    stats.normalizedIntervalsToFindForErrorRate[errorRate]++;
}

// ----------------------------------------------------------------------------
// Function updateMapperStats()
// ----------------------------------------------------------------------------
// Update stats after reading one mapper record.

void updateMapperStats(RabemaOracleStats & stats, Options const & options,
                       BamAlignmentRecord const & mapperRecord, BamAlignmentRecord const & oracleRecord)
{
    SEQAN_ASSERT_EQ(oracleRecord.qName, mapperRecord.qName);

    // Compare mapper record to oracle records.
    if (!hasFlagUnmapped(mapperRecord) && !hasFlagSecondary(mapperRecord))
    {
        unsigned delta = options.delta ? getEditDistance(oracleRecord) : 0;
        unsigned errorRate = getErrorRate(oracleRecord);

        // Update mapper stats.
        stats.mappedReads++;
        stats.readsMappedPerErrorRate[errorRate]++;

        if (isEqual(mapperRecord, oracleRecord, delta))
        {
            // Update mapper stats.
            stats.intervalsFound++;
            stats.normalizedIntervals++;
            stats.intervalsFoundForErrorRate[errorRate]++;
            stats.normalizedIntervalsFoundForErrorRate[errorRate]++;
        }
        else if (options.verbose)
        {
            std::cout << "WARNING: " << mapperRecord.qName
                      << " at " << Pair<unsigned>(mapperRecord.beginPos, getEndPos(mapperRecord))
                      << " differs from " << Pair<unsigned>(oracleRecord.beginPos, getEndPos(oracleRecord)) << '\n';
        }
    }
}

// ----------------------------------------------------------------------------
// Function runRabemaOracle()
// ----------------------------------------------------------------------------

void runRabemaOracle(Options & options)
{
    BamFileIn oracleFile;
    BamFileIn mapperFile;

    BamHeader oracleHeader;
    BamHeader mapperHeader;

    BamAlignmentRecord oracleRecord;
    BamAlignmentRecord mapperRecord;

    RabemaOracleStats stats;

    if (!open(oracleFile, toCString(options.inputOracle)))
        throw RuntimeError("Error while opening the oracle SAM/BAM file.");
    readRecord(oracleHeader, oracleFile);

    if (!open(mapperFile, toCString(options.inputMapper)))
        throw RuntimeError("Error while opening the mapper SAM/BAM file.");
    readRecord(mapperHeader, mapperFile);

    // Use oracle header to interpret the mapper file.
    context(mapperFile) = context(oracleFile);

    // Read first oracle record.
    readRecord(oracleRecord, oracleFile);
    updateOracleStats(stats, oracleRecord);

    // Read first mapper record.
    readRecord(mapperRecord, mapperFile);

    // Read the rest of the oracle and mapper files.
    // Invariant: the oracle file contains each read exactly once.
    do
    {
        bool oracleLTmapper = lessThanSamtoolsQueryName(oracleRecord.qName, mapperRecord.qName);
        bool mapperLToracle = lessThanSamtoolsQueryName(mapperRecord.qName, oracleRecord.qName);

        if (!mapperLToracle && !oracleLTmapper)
            updateMapperStats(stats, options, mapperRecord, oracleRecord);

        if (atEnd(mapperFile) && atEnd(oracleFile))
            break;

        if (!atEnd(oracleFile) && (oracleLTmapper || atEnd(mapperFile)))
        {
            // Read next oracle record.
            readRecord(oracleRecord, oracleFile);
            updateOracleStats(stats, oracleRecord);
        }
        else if (!atEnd(mapperFile))
        {
            // Read next mapper record.
            readRecord(mapperRecord, mapperFile);
        }
    }
    while (true);

    // As precision=found/mapped, just overwrite intervals to find with mapped reads.
    if (options.precision)
    {
        stats.intervalsToFindForErrorRate = stats.readsMappedPerErrorRate;
        stats.normalizedIntervalsToFindForErrorRate = stats.readsMappedPerErrorRate;
    }

    // Write stats to standard out.
    writeStats(std::cout, static_cast<RabemaStats>(stats));

    // Write stats to TSV file.
    if (!empty(options.outputTsv))
    {
        std::ofstream tsvFile(toCString(options.outputTsv), std::ios::out | std::ios::binary);

        if (!tsvFile.good())
            throw RuntimeError("Error while opening the output TSV file.");

        if (writeTsv(tsvFile, static_cast<RabemaStats>(stats), 0u, "oracle", true, "none") != 0)
               throw RuntimeError("Error while writing the output TSV file.");
    }
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    try
    {
        runRabemaOracle(options);
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
