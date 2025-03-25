//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include<bit>
#include "util.hpp"
#include "../KMC/kmc_api/kmer_api.h"
#include "../KMC/kmc_api/kmc_file.h"
#include "htslib/bgzf.h"
//#include "nc_utils.h"
//
//
#define BIN_LEN 200000

struct KMCdb {
    std::string root;
    std::vector<CKMCFile*> files;
    int ngenomes, nbytes, kmer_length;

    KMCdb(std::string root_, int ngenomes_) : root(root_), ngenomes(ngenomes_) {
        size_t ndbs = ceil(((float)ngenomes) / 32);

        std::stringstream ss;
        files.resize(ndbs);
        for (size_t i = 0; i < ndbs; i++) {
            ss << root << "/kmc/bitvec" << i;
            files[i] = new CKMCFile();
            files[i]->OpenForRA(ss.str());
            ss.str("");
            kmer_length = files[i]->KmerLength();
        }

        nbytes = (int) ceil(((float)ngenomes ) / 8);
    }

    void anchor_fasta(std::string name, std::string fasta) {
        std::ifstream input_file;
        input_file.open(fasta);

        std::stringstream ss;

        ss << root << "/anchor/" << name << "/bitmap.1";
        auto basename = ss.str();
        ss << ".gz";
        auto bitmap1 = bgzf_open(ss.str().c_str(), "w");
        bgzf_index_build_init(bitmap1);
        ss.str("");

        ss << root << "/anchor/" << name << "/bitmap.100";
        auto basename100 = ss.str();
        ss << ".gz";
        auto bitmap100 = bgzf_open(ss.str().c_str(), "w");
        bgzf_index_build_init(bitmap100);
        ss.str("");

        ss << root << "/anchor/" << name << "/bitsum.bins.tsv";
        auto bitsum_out = std::ofstream(ss.str());
        bitsum_out << "chr\tstart";
        for (size_t i = 0; i < ngenomes+1; i++) {
            bitsum_out << "\t" << i;
        }
        bitsum_out << "\n";
        ss.str("");

        ss << root << "/anchor/" << name << "/chrs.tsv";
        auto chrs_out = std::ofstream(ss.str());
        ss.str("");
        chrs_out << "name\tid\tsize\tgene_count\n";

        size_t offset = 0;
        size_t chr_num = 0;

        if (input_file.is_open()) {
            std::string line, name, seq;

            while (std::getline(input_file, line)) {
                if (line[0] == '>') {
                    if (!name.empty()) {
                        seq = ss.str();
                        auto ksize = write_bits(chr_num, seq, bitmap1, bitmap100, bitsum_out);
                        ss.str("");

                        size_t s = name.find(' ');
                        chrs_out << name.substr(0,s) << "\t" << chr_num << "\t" << ksize << "\t0\n";
                        chrs_out.flush();
                        chr_num += 1;
                    }
                    name = line.substr(1);
                } else {
                    ss << line;
                }
            }

            seq = ss.str();
            auto ksize = write_bits(chr_num, seq, bitmap1, bitmap100, bitsum_out);
            size_t s = name.find(' ');
            chrs_out << name.substr(0,s) << "\t" << chr_num << "\t" << ksize << "\t0\n";
            chrs_out.flush();
        }

        bgzf_index_dump(bitmap1, basename.c_str(), ".gzi");
        bgzf_close(bitmap1);

        bgzf_index_dump(bitmap100, basename100.c_str(), ".gzi");
        bgzf_close(bitmap100);
        bitsum_out.close();
        chrs_out.close();
    }

    //void write_bits(std::string &seq, std::ofstream &bitmap1) {
    size_t write_bits(size_t chr_num, std::string &seq, BGZF *bitmap1, BGZF *bitmap100, std::ofstream &bitsum_out) {

        size_t binlen = 200000;
        size_t nkmers = seq.size() - kmer_length + 1;
        if (nkmers / binlen < 100) {
            binlen = nkmers / 100;
        }

        i64 nchunks = (nkmers / binlen) + (nkmers % binlen != 0);
        i64 chunk_start = 0;

        for (size_t chunk = 0; chunk < nchunks; chunk++) {
            auto chunk_end = std::min(chunk_start+binlen, nkmers);


            auto seq_chunk = seq.substr(chunk_start, (chunk_end-chunk_start) + kmer_length - 1);
            //std::cout << chunk_start << " " << chunk_end << " " << seq_chunk.size() << "\n";


            size_t dbi = 0;
            size_t n = nbytes;

            std::vector<u8> bytes;
            std::vector<int> popcnts;

            size_t offs = 0;
            for (size_t dbi = 0; dbi < files.size(); dbi++) {
                if (nbytes <= 4) {
                    n = nbytes;
                } else if (dbi == files.size()-1 and nbytes % 4 > 0) {
                    n = nbytes % 4;
                } else {
                    n = 4;
                }

                std::vector<u32> ints;
                files[dbi]->GetCountersForRead(seq_chunk, ints);

                if (bytes.size() == 0) {
                    bytes.resize(ints.size()*nbytes);
                    popcnts.resize(ints.size());
                }
                
                size_t i = offs;
                for (size_t j = 0; j < ints.size(); j++) {
                    for (size_t sh = 0; sh < n; sh++) {
                        bytes[i+sh] = (u8) (ints[j] >> (sh*8));
                    }
                    popcnts[j] += __builtin_popcount(ints[j]);
                    i += nbytes;
                }

                offs += n;
            }

            bgzf_write(bitmap1, (void*) bytes.data(), bytes.size() * sizeof(u8));

            size_t i100 = 0; //ceil( ((float) ((bytes.size() / nbytes)) / 100 ));
            for (size_t i = ((100-(chunk_start % 100))%100)*nbytes; i < bytes.size(); i += 100*nbytes) {
                for (size_t b = 0; b < nbytes; b++) {
                    bytes[i100] = bytes[i+b];
                    i100 += 1;
                }
            }
            bytes.resize(i100);
            bgzf_write(bitmap100, (void*) bytes.data(), bytes.size() * sizeof(u8));

            std::vector<size_t> bitsums;
            bitsums.resize(ngenomes+1);
            for (size_t k = 0; k < binlen && k < popcnts.size(); k++) {
                bitsums[popcnts[k]]++;
            }
            bitsum_out << chr_num << "\t" << chunk_start;
            for (auto c : bitsums) {
                bitsum_out << "\t" << c;
            }
            bitsum_out << "\n";
            bitsum_out.flush();
            chunk_start += binlen;
        }

        return nkmers;

    }

    ~KMCdb() {
        for (auto p : files) {
            delete p;
        }
    }
};

int main(int argc, char* argv[]) {
	int ngenomes = atoi(argv[1]);
	std::string root(argv[2]); //, fasta(argv[3]), outfile(argv[4]);

    if (argc > (ngenomes*2)+3) {
        std::cerr << "Error: expected " << (ngenomes*2+1) << " or fewer arguments\n";
        return 1;
    }

    KMCdb db(root, ngenomes);

    #pragma omp parallel for
    for (size_t i = 3; i < argc; i+=2) {
        auto name = std::string(argv[i]), 
             fasta = std::string(argv[i+1]);
        std::cout << "Anchoring " << name << " " << fasta << "\n";
        db.anchor_fasta(name, fasta);
    }

    //currently assuming that there is only a single sequence
    //and that it is on only one line of the file

	//32 indicates bits per k-mer, which could be more compressed in the future


	//ann_out.close();
    return 0;
}

