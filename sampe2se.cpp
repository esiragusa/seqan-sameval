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

#include <seqan/bam_io.h>

using namespace seqan;

// ============================================================================
// Functors
// ============================================================================

typedef EqualsChar<'.'>        IsDot;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clearPairedFlags()
// ----------------------------------------------------------------------------

void clearPairedFlags(BamAlignmentRecord & record)
{
    if (hasFlagMultiple(record))
        record.flag = record.flag ^ BAM_FLAG_MULTIPLE;
    if (hasFlagFirst(record))
        record.flag = record.flag ^ BAM_FLAG_FIRST;
    if (hasFlagLast(record))
        record.flag = record.flag ^ BAM_FLAG_LAST;
    if (hasFlagNextRC(record))
        record.flag = record.flag ^ BAM_FLAG_NEXT_RC;
    if (hasFlagNextUnmapped(record))
        record.flag = record.flag ^ BAM_FLAG_NEXT_UNMAPPED;
    if (hasFlagAllProper(record))
        record.flag = record.flag ^ BAM_FLAG_ALL_PROPER;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cerr << "usage: sampe2se <pe-bam> <se-bam>" << std::endl;
        return 1;
    }

    BamFileIn bamFileIn;
    BamFileOut bamFileOut(bamFileIn);
    BamHeader header;
    BamAlignmentRecord record;

    if (!open(bamFileIn, toCString(argv[1])))
    {
        std::cerr << "Can't open the input file." << std::endl;
        return 1;
    }

    if (!open(bamFileOut, toCString(argv[2])))
    {
        std::cerr << "Can't open the output file." << std::endl;
        return 1;
    }

    readRecord(header, bamFileIn);
    writeRecord(bamFileOut, header);

    StringSet<CharString> qName;

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);

        if (hasFlagFirst(record) || endsWith(record.qName, "/1"))
        {
            if (endsWith(record.qName, "/1"))
                resize(record.qName, length(record.qName) - 2);

            clear(qName);
            strSplit(qName, record.qName, IsDot());
            clear(record.qName);
            append(record.qName, qName[0]);
            appendValue(record.qName, '.');
            appendValue(record.qName, 'L');
            append(record.qName, qName[1]);
        }
        else if (hasFlagLast(record) || endsWith(record.qName, "/2"))
        {
            if (endsWith(record.qName, "/2"))
                resize(record.qName, length(record.qName) - 2);

            clear(qName);
            strSplit(qName, record.qName, IsDot());
            clear(record.qName);
            append(record.qName, qName[0]);
            appendValue(record.qName, '.');
            appendValue(record.qName, 'R');
            append(record.qName, qName[1]);
        }

        clearPairedFlags(record);

        writeRecord(bamFileOut, record);
    }

    return 0;
}

