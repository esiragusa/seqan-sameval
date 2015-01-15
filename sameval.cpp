// ==========================================================================
//                               SamEval
// ==========================================================================
// Copyright (c) 2011-2015, Enrico Siragusa, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

using namespace seqan;

// ============================================================================
// Options
// ============================================================================

struct Options
{
    CharString inputOracle;
    CharString inputMapper;
    CharString outputTsv;

    bool mappingQuality;
    bool verbose;
    bool delta;

    Options() :
        mappingQuality(false),
        verbose(false),
        delta(false)
    {}
};

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------

void setupArgumentParser(ArgumentParser & parser)
{
    setAppName(parser, "sameval");
    setShortDescription(parser, "SAMeval");
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
    addOption(parser, seqan::ArgParseOption("q", "mapping-quality", "Class by mapping quality (QUAL). \
                                                                     Default: class by oracle error rate (NM)."));
    addOption(parser, seqan::ArgParseOption("d", "delta", "Tolerate +- delta bp for begin/end positions. \
                                                           Default: require exact begin/end positions."));

    addSection(parser, "Input / Output");

    addOption(parser, seqan::ArgParseOption("", "out-tsv", "Write the statistics to this TSV file.",
                                            seqan::ArgParseArgument::OUTPUT_FILE, "TSV"));
    setValidValues(parser, "out-tsv", "tsv");
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
    getOptionValue(options.mappingQuality, parser, "mapping-quality");
    getOptionValue(options.delta, parser, "delta");
    getOptionValue(options.outputTsv, parser, "out-tsv");

    return seqan::ArgumentParser::PARSE_OK;
}


// ============================================================================
// Functions
// ============================================================================

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
    SEQAN_ASSERT_GEQ(value, (TValue)delta);
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

// ============================================================================
// Stats
// ============================================================================

// ----------------------------------------------------------------------------
// Class Stats
// ----------------------------------------------------------------------------

typedef Triple<__uint64>    Counts;
typedef String<Counts>      Stats;

// ----------------------------------------------------------------------------
// Accessors
// ----------------------------------------------------------------------------

template <typename T1, typename T2, typename T3, typename TSpec>
inline T1 &
totalCount(Triple<T1, T2, T3, TSpec> & triple)
{
    return triple.i1;
}

template <typename T1, typename T2, typename T3, typename TSpec>
inline T1 const &
totalCount(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i1;
}

template <typename T1, typename T2, typename T3, typename TSpec>
inline T2 &
correctCount(Triple<T1, T2, T3, TSpec> & triple)
{
    return triple.i2;
}

template <typename T1, typename T2, typename T3, typename TSpec>
inline T2 const &
correctCount(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i2;
}

template <typename T1, typename T2, typename T3, typename TSpec>
inline T3 &
incorrectCount(Triple<T1, T2, T3, TSpec> & triple)
{
    return triple.i3;
}

template <typename T1, typename T2, typename T3, typename TSpec>
inline T3 const &
incorrectCount(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i3;
}

// ----------------------------------------------------------------------------
// Function resizeStats()
// ----------------------------------------------------------------------------

void resizeStats(Stats & stats, unsigned threshold)
{
    if (length(stats) <= threshold + 1)
        resize(stats, threshold + 1, Counts(0, 0, 0));
}

// ----------------------------------------------------------------------------
// Function writeStats()
// ----------------------------------------------------------------------------

template <typename TStream>
void writeStats(TStream & stream, Stats const & stats)
{
    stream << "threshold" << '\t' << "total" << '\t' << "correct" << '\t' << "incorrect" << '\n';

    for (unsigned i = 0; i < length(stats); ++i)
        stream << i << '\t' << totalCount(stats[i]) << '\t' << correctCount(stats[i]) << '\t' << incorrectCount(stats[i]) << '\n';
}

// ----------------------------------------------------------------------------
// Function updateOracleStats()
// ----------------------------------------------------------------------------
// Update stats after reading one oracle record.

void updateOracleStats(Stats & stats, Options const & options, BamAlignmentRecord const & record)
{
    if (!options.mappingQuality)
    {
        unsigned errorRate = getErrorRate(record);
        resizeStats(stats, errorRate);
        totalCount(stats[errorRate])++;
    }
}

// ----------------------------------------------------------------------------
// Function updateMapperStats()
// ----------------------------------------------------------------------------
// Update stats after reading one mapper record.

void updateMapperStats(Stats & stats, Options const & options, BamAlignmentRecord const & record)
{
    if (options.mappingQuality)
    {
        unsigned mappingQuality = record.mapQ;
        resizeStats(stats, mappingQuality);
        totalCount(stats[mappingQuality])++;
    }
}

// ----------------------------------------------------------------------------
// Function updateStats()
// ----------------------------------------------------------------------------
// Update stats after reading one mapper and oracle record.

void updateStats(Stats & stats, Options const & options,
                 BamAlignmentRecord const & mapperRecord, BamAlignmentRecord const & oracleRecord)
{
    SEQAN_ASSERT_EQ(oracleRecord.qName, mapperRecord.qName);

    // Compare mapper record to oracle records.
    if (!hasFlagUnmapped(mapperRecord) && !hasFlagSecondary(mapperRecord))
    {
        unsigned delta = options.delta ? getEditDistance(oracleRecord) : 0;
        unsigned errorRate = getErrorRate(oracleRecord);
        unsigned threshold = options.mappingQuality ? mapperRecord.mapQ : errorRate;

        // Update mapper stats.
        if (isEqual(mapperRecord, oracleRecord, delta))
        {
            correctCount(stats[threshold])++;
        }
        else
        {
            incorrectCount(stats[threshold])++;

            if (options.verbose)
            {
                std::cerr << "WARNING: " << mapperRecord.qName
                          << " at " << Triple<unsigned>(mapperRecord.beginPos, getEndPos(mapperRecord), getEditDistance(mapperRecord))
                          << " instead of " << Triple<unsigned>(oracleRecord.beginPos, getEndPos(oracleRecord), getEditDistance(oracleRecord)) << '\n';
            }
        }
    }
}

// ============================================================================
// Main
// ============================================================================

// ----------------------------------------------------------------------------
// Function run()
// ----------------------------------------------------------------------------

void run(Options & options)
{
    BamFileIn oracleFile;
    BamFileIn mapperFile;

    BamHeader oracleHeader;
    BamHeader mapperHeader;

    BamAlignmentRecord oracleRecord;
    BamAlignmentRecord mapperRecord;

    Stats stats;

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
    updateOracleStats(stats, options, oracleRecord);

    // Read first mapper record.
    readRecord(mapperRecord, mapperFile);
    updateMapperStats(stats, options, mapperRecord);

    // Read the rest of the oracle and mapper files.
    // Invariant: the oracle file contains each read exactly once.
    do
    {
        bool oracleLTmapper = lessThanSamtoolsQueryName(oracleRecord.qName, mapperRecord.qName);
        bool mapperLToracle = lessThanSamtoolsQueryName(mapperRecord.qName, oracleRecord.qName);

        if (!mapperLToracle && !oracleLTmapper)
            updateStats(stats, options, mapperRecord, oracleRecord);

        if (atEnd(mapperFile) && atEnd(oracleFile))
            break;

        if (!atEnd(oracleFile) && (oracleLTmapper || atEnd(mapperFile)))
        {
            // Read next oracle record.
            readRecord(oracleRecord, oracleFile);
            updateOracleStats(stats, options, oracleRecord);
        }
        else if (!atEnd(mapperFile))
        {
            // Read next mapper record.
            readRecord(mapperRecord, mapperFile);
            updateMapperStats(stats, options, mapperRecord);
        }
    }
    while (true);

    // Write stats to standard out.
    writeStats(std::cout, stats);

    // Write stats to TSV file.
    if (!empty(options.outputTsv))
    {
        std::ofstream tsvFile;

        if (!open(tsvFile, toCString(options.outputTsv)))
            throw RuntimeError("Error while opening the output TSV file.");

        writeStats(tsvFile, stats);
    }
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser);

    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    try
    {
        run(options);
    }
    catch (Exception const & e)
    {
        std::cerr << getAppName(parser) << ": " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
