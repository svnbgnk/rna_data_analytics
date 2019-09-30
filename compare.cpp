#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>
#include <seqan/bam_io.h>
#include<iostream>
#include<fstream>
#include <seqan/find.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <map>
#include <string>

using namespace seqan;
using namespace std;

//TODO star did not write down unmapped reads to to rerun script

struct Stats
{
    uint32_t reads_count = 0;
    uint32_t missing1 = 0;
    uint32_t missing2 = 0;


    uint32_t ambiguous_count = 0;
    uint32_t ambiguous1_count = 0;
    uint32_t ambiguous2_count = 0;

    uint32_t unmapped_count = 0;
    uint32_t unmapped1_count = 0;
    uint32_t unmapped2_count = 0;
    uint32_t similar_count = 0;
    uint32_t unsimilar_count = 0;
    std::vector<std::pair<CharString, uint32_t> > differences;
    std::vector<std::pair<CharString, int> > eventsdif;

    void print()
    {
        std::cerr << "Overall reads: " << 2*reads_count << "\n";

        std::cerr << "Global ambiguous: " << ambiguous_count << "\n";
        std::cerr << "Global unmapped: " << unmapped_count << "\n";

        std::cerr << "1 ambiguous: " << ambiguous1_count << "\n";
        std::cerr << "1 unmapped: " << unmapped1_count << "\n";

        std::cerr << "2 ambiguous: " << ambiguous2_count << "\n";
        std::cerr << "2 unmapped: " << unmapped2_count << "\n";

        std::cerr << "Similar: " << similar_count << "\n";
    }

    void sim_diff()
    {
        std::cerr << "Cigar of similar alignments differ by: \n";
        for(int i = 0; i < differences.size(); ++i){
           std::cerr << "<" << toCString(differences[i].first) << ", " << differences[i].second << ">, ";
           if(i % 10 == 0 && i != 0)
               std::cerr << "\n";
        }
    }

    void sim_events()
    {
        std::cerr << "Events of similar alignments differ by: \n";
        for(int i = 0; i < eventsdif.size(); ++i){
           std::cerr << "<" << toCString(eventsdif[i].first) << ", " << eventsdif[i].second << ">, ";
           if(i % 10 == 0 && i != 0)
               std::cerr << "\n";
        }
    }
};


void load_alignments_of_read(BamFileIn & bamFileIn1,
                             BamFileIn & bamFileIn2,
                             auto & records1,
                             auto & records2,
                             BamAlignmentRecord & lastRecord1,
                             BamAlignmentRecord & lastRecord2,
                             auto & stats)
{

    BamAlignmentRecord record1;
    BamAlignmentRecord record2;
    std::string l1 = toCString(lastRecord1.qName);
    std::string l2 = toCString(lastRecord2.qName);

    while(l1.compare(l2) != 0){

        if(l1.compare(l2) < 0){
            readRecord(lastRecord1, bamFileIn1);
            l1 = toCString(lastRecord1.qName);
            ++stats.missing1;
            std::cerr << "Missing in second file: " << l1 << "\n";
        }
        else
        {
            readRecord(lastRecord2, bamFileIn2);
            l2 = toCString(lastRecord2.qName);
            std::cerr << "Missing in first file: " << l2 << "\n";
            ++stats.missing2;
        }
    }

    std::cerr << "r1: " << toCString(lastRecord1.qName) << "\n";
    std::cerr << "r2: " << toCString(lastRecord2.qName) << "\n";

    records1.push_back(lastRecord1);
    records2.push_back(lastRecord2);

    while (!atEnd(bamFileIn1))
    {
        try
        {
            readRecord(lastRecord1, bamFileIn1);
            std::string l1_tmp = toCString(lastRecord1.qName);
            if(l1.compare(l1_tmp) != 0)
                break;
            records1.push_back(lastRecord1);

        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }

    while (!atEnd(bamFileIn2))
    {
        try
        {
            readRecord(lastRecord2, bamFileIn2);
            std::string l2_tmp = toCString(lastRecord2.qName);
            if(l2.compare(l2_tmp) != 0)
                break;
            records2.push_back(lastRecord2);

        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }
    if(atEnd(bamFileIn1)){
        lastRecord1.qName = "zend0000";
    }

    if(atEnd(bamFileIn2)){
        lastRecord2.qName = "zend0000";
    }
}

template<typename T>
void printv (std::vector <T> v)
{
    std::cerr << toCString(v[0].qName) << "\n";
    for(int i = 0; i < v.size(); ++i){
        std::cerr << v[i].beginPos << ", ";
    }
    std::cerr << "\n";
}

template<typename TStream, typename TCigarString>
void print_cigar(TStream & stream, TCigarString const cigar)
{
    for(int i = 0; i < length(cigar); ++i){
        stream << (int)cigar[i].count << (char)cigar[i].operation;
    }
}

template <typename TStream, typename TRecord>
inline void printRecord(TStream & stream, TRecord const & me)
{
    CharString seq = me.seq;
    CharString tags = me.tags;
    stream << toCString(me.qName) << "\t" << me.flag << "\t" << me.rID << "\t" << me.beginPos << "\t" << me.mapQ << "\t" << me.bin << "\t";
    print_cigar(stream, me.cigar);
    stream << "\t" << me.rNextId << "\t" << me.pNext << "\t" << me.tLen  << "\t" << toCString(seq)  << "\t" << toCString(me.qual) << "\t" << toCString(tags)<< "\n";
}


bool posSmaller(BamAlignmentRecord & x, BamAlignmentRecord const & y)
{
    return (x.beginPos < y.beginPos);
}

template <typename TStream, typename TRecord>
inline void printRecords(TStream & stream, std::vector<TRecord> const & records)
{
    for(int i = 0; i < records.size(); ++i){
        printRecord(stream, records[i]);
    }
}

int softThreshold = 5;
int matchThreshold = 15;
int insertionThreshold = 15;
CigarElement<> soft = seqan::CigarElement<>('S', 1);
CigarElement<> match = seqan::CigarElement<>('M', 1);
CigarElement<> del = seqan::CigarElement<>('D', 1);
CigarElement<> ins = seqan::CigarElement<>('I', 1);

int detectevents(auto & cigar1){
    int event = 0;
    int prefixMatches = 0;
    int op = 0;
    if(soft.operation == cigar1[op].operation && cigar1[op].count > softThreshold)
    {
        ++event;
    }
    for(; op < length(cigar1); ++op){
        if(match.operation == cigar1[op].operation)
        {
            prefixMatches += cigar1[op].count;
        }
        else if(prefixMatches > matchThreshold && ins.operation == cigar1[op].operation && insertionThreshold < cigar1[op].count)
        {
            ++event;
            prefixMatches = 0;
        }
    }
    --op;

    if(soft.operation == cigar1[op].operation && cigar1[op].count > softThreshold)
    {
        ++event;
    }
    return event;
}


int main(int argc, char const * argv[])
{
    ArgumentParser parser("Search");

    addOption(parser, ArgParseOption("s1", "sam1", "Path to the sam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "sam1");

    addOption(parser, ArgParseOption("s2", "sam2", "Path to the sam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "sam2");

//     setRequired(parser, "reads1");

    addOption(parser, ArgParseOption("o", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
//     setRequired(parser, "output");

    addOption(parser, ArgParseOption("bS", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString samPath1, samPath2, outputPath;
    int batchSize = 100000, barcodeLength;

    getOptionValue(samPath1, parser, "sam1");
    getOptionValue(samPath2, parser, "sam2");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(batchSize, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");

    //prepare loading sam
    BamFileIn bamFileIn1;
    BamHeader header1;
    if (!open(bamFileIn1, toCString(samPath1)))
    {
        std::cerr << "ERROR: could not open input file " << samPath1 << ".\n";
        return 1;
    }
    readHeader(header1, bamFileIn1);
//     std::cerr << "check2\n";

    //prepare loading second sam
    BamFileIn bamFileIn2;
    BamHeader header2;
    if (!open(bamFileIn2, toCString(samPath2)))
    {
        std::cerr << "ERROR: could not open input file " << samPath2 << ".\n";
        return 1;
    }
    readHeader(header2, bamFileIn2);

    Stats stats;

    bool eventdetection = true;


    BamAlignmentRecord lastRecord1;
    BamAlignmentRecord lastRecord2;
    lastRecord1.qName = "00NA1";
    lastRecord2.qName = "00NA2";

    std::string l1, l2;

    while(l1.compare("zend0000") != 0 && l2.compare("zend0000") != 0)
    {
        ++stats.reads_count;
        l1 = toCString(lastRecord1.qName);
        l2 = toCString(lastRecord2.qName);
        std::vector<BamAlignmentRecord> records1;
        std::vector<BamAlignmentRecord> records2;


        load_alignments_of_read(bamFileIn1, bamFileIn2, records1, records2, lastRecord1, lastRecord2, stats);

        std::cerr << "Read " << toCString(records1[0].qName) << "\t" << toCString(records2[0].qName) << "\t" << stats.reads_count <<  "-------------------\n";
        std::cerr << "records size: " << records1.size() << "\t" << records2.size() << "\n";

//         if(stats.reads_count > 20)
//             exit(0);

        bool single = false;

        if(records1.size() < 2 || records2.size() < 2)
        {
           std::cerr << "Single Read??\n";
            if(records1.size() == 1 && records2.size() == 1)
                single = true;
            printv(records1);
            std::cerr << "\n";
            printv(records2);
        }

        if(records1.size() > 2 || records2.size() > 2)
        {
            std::sort(records1.begin(), records1.end(), posSmaller);
            std::sort(records2.begin(), records2.end(), posSmaller);
            std::cerr << "Multimappers\n";
            std::cerr << "File 1: \n";
            printRecords(std::cerr, records1);
            std::cerr << "File 2: \n";
            printRecords(std::cerr, records2);

            if(hasFlagMultiple(records1[0]) || hasFlagMultiple(records1[1]))
            {
                ++stats.ambiguous1_count;
            }
            if(hasFlagMultiple(records2[0]) || hasFlagMultiple(records2[1]))
            {
                ++stats.ambiguous2_count;
            }

            if(records1.size() == 2 || records2.size() == 2)
            {
                bool first = records1.size() == 2;
                std::vector<BamAlignmentRecord> & single = (first) ? records1 : records2;
                std::vector<BamAlignmentRecord> & multi = (!first) ? records1 : records2;
                std::vector<BamAlignmentRecord> newMulti;
                for(int j = 0; j < 2; ++j){
                    for(int i = 0; i < multi.size(); ++i){
                        bool similar = (single[j].beginPos < (multi[i].beginPos + 10) && single[j].beginPos + 10 > (multi[i].beginPos));
                        std::cout << "Check: " << single[j].beginPos << "\t" << multi[i].beginPos << "\n";
                        if(similar){
                            std::cout << "Sim\n";
                            newMulti.push_back(multi[i]);
                        }
                    }
                }
                std::cout << "multi size: " << newMulti.size() << "\n";
                if (newMulti.size() == 2)
                {
                    multi = newMulti;
                    std::cerr << "Size " << multi.size() << "\t" << records1.size() << "\t" << records2.size() << "\n";
                    std::cerr << "rescued unique alignment from multimapper\n";
                }
            }
        }


        if(records1.size() == 2 && records2.size() == 2)
        {
            std::cerr << "Paired uniquely mapping reads\n";

            if(hasFlagUnmapped(records1[0]) || hasFlagUnmapped(records1[1]))
            {
                ++stats.unmapped1_count;

            }

            if(hasFlagUnmapped(records2[0]) || hasFlagUnmapped(records2[1]))
            {
                ++stats.unmapped2_count;
            }

            if(hasFlagUnmapped(records1[0]) || hasFlagUnmapped(records1[1]) || hasFlagUnmapped(records2[0]) || hasFlagUnmapped(records2[1]))
            {
                ++stats.unmapped_count;
            }
            else
            {

    //             std::cerr << "Order\n";

                if(hasFlagFirst(records1[1])){
                    std::swap(records1[0], records1[1]);
                }
                if(hasFlagFirst(records2[1])){
                    std::swap(records2[0], records2[1]);
                }
    //             printv(records1);
    //             std::cerr << "\n";
    //             printv(records2);


                for(int i = 0; i < (2 - single); ++i){
                    bool similar = (records1[i].beginPos < (records2[i].beginPos + 10) && records1[i].beginPos + 10 > (records2[i].beginPos));
                    stats.similar_count += similar;

                    if(similar)
                    {
    //                     std::cerr << "Check Cigar\n";
                        String<CigarElement<> > cigar1 = records1[i].cigar;
                        String<CigarElement<> > cigar2 = records2[i].cigar;
                        int16_t pos1 = 0, pos2 = 0, opos1 = 0, opos2 = 0;

                        uint32_t dif = 0;
                        uint32_t cigarpos = 0;
                        while(true)
                        {

                            if(pos1 >= length(cigar1) && pos2 >= length(cigar2))
                            {
                                break;
                            }
                            else if(pos1 >= length(cigar1) || pos2 >= length(cigar2))
                            {
                                ++dif;
                            }
                            else if (cigar1[pos1].operation != cigar2[pos2].operation)
                            {
                                ++dif;
                            }
    //                         std::cerr << cigarpos << "\t" << dif << "\t" << cigar1[pos1].operation << "\t" << cigar2[pos2].operation << "\n";
                            //TODO jump muliple position in a single operation
                            ++opos1;
                            ++opos2;
                            if(opos1 >= cigar1[pos1].count){
                                ++pos1;
                                opos1 = 0;
                            }

                            if(opos2 >= cigar2[pos2].count){
                                ++pos2;
                                opos2 = 0;
                            }

    //                         std::cerr << "counts " << cigar1[pos1].count << "\t" << cigar2[pos2].count

                            ++cigarpos;
                        }

                        stats.differences.push_back(std::make_pair(toCString(records1[i].qName), dif));
                        if(dif == 0){
    //                         std::cerr << "same Cigar\n";
                        }
                        else
                        {
                            std::cerr << "Diff: " << dif << "\n";
                            printRecord(std::cerr, records1[i]);
                            printRecord(std::cerr, records2[i]);
                        }


                        if(eventdetection)
                        {
                            int events = detectevents(cigar1);
                            int events2 = detectevents(cigar2);
                            stats.eventsdif.push_back(std::make_pair(toCString(records1[i].qName), (events - events2)));

                        }
                    }
                    else
                    {
                        ++stats.unsimilar_count;
                    }

                }
            }

            //TODO Compare eventdetection
        }
    }



    std::cerr << "\n";
    stats.print();
    stats.sim_diff();

    std::cerr << "Finished!\n";

    return 0;
}



