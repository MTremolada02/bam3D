#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <map>
#include <algorithm>

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include "functions.h"
#include "global.h"
#include "def.hpp"

void Runner::loadInput(UserInputBam3D userInput) {
    this->userInput = userInput;
}

void Runner::write_section_header(std::ofstream& myfile, const std::string& section_name, const std::string& columns_line)
{
    myfile << "\n#" << section_name << "\n";
    myfile << columns_line << "\n";
}

void Runner::write_all_stats_file(const std::string& out_path,bam_hdr_t* bamHdr)
{
    std::ofstream myfile(out_path, std::ios::out);

    if (!myfile.is_open()) {
        std::cout << "cannot open " << out_path << std::endl;
        return;
    }

	if(userInput.pair_read_stats){ 

		if (!graph.graph_inter_count.empty() && !graph.graph_intra_count.empty()) {
			write_reference_graph(myfile,bamHdr);
		}

		if (!graph.Ps_binned_dist_count.empty()) {
			write_binned_map(myfile, "DIST_PS", graph.Ps_binned_dist_count);
		}

		if (!graph.Ps_binned_dist_count.empty()) {
			write_binned_map(myfile, "DIST_FF", graph.ff_binned_dist_count);
		}

		if (!graph.Ps_binned_dist_count.empty()) {
			write_binned_map(myfile, "DIST_FR", graph.fr_binned_dist_count);
		}

		if (!graph.Ps_binned_dist_count.empty()) {
			write_binned_map(myfile, "DIST_RF", graph.rf_binned_dist_count);
		}


		if (!graph.Ps_binned_dist_count.empty()) {
			write_binned_map(myfile, "DIST_RR", graph.rr_binned_dist_count);
		}

		write_pair_types_section(myfile);
	}

    myfile.close();
}


//GRAPHS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Runner::write_reference_graph(std::ofstream& myfile, bam_hdr_t* bamHdr) {
    //nodes
    write_section_header(myfile, "GRAPH_NODES", "id\tname\tlength\tintra_count");

    for (int32_t tid = 0; tid < bamHdr->n_targets; ++tid) {
        uint64_t intra = 0;
        auto it = graph.graph_intra_count.find(tid);
        if (it != graph.graph_intra_count.end()) intra = it->second;

        myfile << tid << "\t" << bamHdr->target_name[tid] << "\t" << bamHdr->target_len[tid] << "\t"<< intra << "\n";
    }
    
    // archi
    write_section_header(myfile, "GRAPH_EDGES", "source\ttarget\tinter_count");

    for (const auto& kv : graph.graph_inter_count) {
        uint32_t source = kv.first.first;
        uint32_t target = kv.first.second;
        uint64_t count  = kv.second;

        myfile << source << "\t" << target << "\t" << count << "\n";
    }
}

void Runner::write_binned_map(std::ofstream& myfile, const std::string& section_name, const std::unordered_map<uint64_t, uint64_t>& data_map)
{
    write_section_header(myfile, section_name, "bin_start\tbin_end\tcount\tcount_per_bp\tcount_fraction");

    std::vector<std::pair<uint32_t, uint64_t>> entries(data_map.begin(), data_map.end());

    std::sort(entries.begin(), entries.end(),
              [](const auto& a, const auto& b) {
                  return a.first < b.first;
              });

	uint64_t total_count = 0;
	for (const auto& kv : entries) {
  		total_count += kv.second;
	}

    for (const auto& kv : entries) {
        uint32_t bin_index = kv.first;
        uint64_t count = kv.second;

        uint64_t bin_start = static_cast<uint64_t>(std::floor(std::pow(graph.log_bin_factor, bin_index)));
        uint64_t bin_end   = static_cast<uint64_t>(std::floor(std::pow(graph.log_bin_factor, bin_index + 1))) - 1;

        if (bin_start < 1) bin_start = 1;
        if (bin_end < bin_start) bin_end = bin_start;

		uint64_t bin_width = bin_end - bin_start + 1;
		double count_per_bp = static_cast<double>(count) / static_cast<double>(bin_width);

		double count_fraction = 0.0;
		if (total_count > 0) {
			count_fraction = static_cast<double>(count) / static_cast<double>(total_count);
		}

        myfile << bin_start << "\t" << bin_end << "\t" << count << "\t" << count_per_bp << "\t" << count_fraction << "\n";
    }
}

void Runner::update_log_binned_distance(uint64_t dist, std::unordered_map<uint64_t, uint64_t>& specific_map)
{
    if (dist == 0) return;

    uint32_t bin_index = static_cast<uint32_t>(
        std::floor(std::log((double)dist) * graph.inv_log_bin_factor) //std::log((double)dist): Trasforma la distanza in una scala logaritmica, / std::log(factor) (o * inv_factor): Serve a decidere la "larghezza" dei tuoi scalini. Se il factor è piccolo (es. 1.01), avrai tantissimi scalini stretti
    );

    ++specific_map[bin_index];
}

void Runner::write_pair_types_section(std::ofstream& myfile)
{
    write_section_header(
        myfile,
        "PAIR_TYPES",
        "run\tUU\tRU\tUR\tWW\tDD\tMU\tMR\tMM\tNM\tNU\tNR\tNN\tpUU\tpRU\tpUR\tpWW\tpDD\tpMU\tpMR\tpMM\tpNM\tpNU\tpNR\tpNN"
    );

    const double total_pairs =
        qnameStats.UU + qnameStats.RU + qnameStats.UR + qnameStats.WW +
        qnameStats.DD + qnameStats.MU + qnameStats.MR + qnameStats.MM +
        qnameStats.NM + qnameStats.NU + qnameStats.NR + qnameStats.NN;

    myfile << "run1" << "\t"
           << qnameStats.UU << "\t"
           << qnameStats.RU << "\t"
           << qnameStats.UR << "\t"
           << qnameStats.WW << "\t"
           << qnameStats.DD << "\t"
           << qnameStats.MU << "\t"
           << qnameStats.MR << "\t"
           << qnameStats.MM << "\t"
           << qnameStats.NM << "\t"
           << qnameStats.NU << "\t"
           << qnameStats.NR << "\t"
           << qnameStats.NN << "\t"
           << percentage(qnameStats.UU, total_pairs) << "\t"
           << percentage(qnameStats.RU, total_pairs) << "\t"
           << percentage(qnameStats.UR, total_pairs) << "\t"
           << percentage(qnameStats.WW, total_pairs) << "\t"
           << percentage(qnameStats.DD, total_pairs) << "\t"
           << percentage(qnameStats.MU, total_pairs) << "\t"
           << percentage(qnameStats.MR, total_pairs) << "\t"
           << percentage(qnameStats.MM, total_pairs) << "\t"
           << percentage(qnameStats.NM, total_pairs) << "\t"
           << percentage(qnameStats.NU, total_pairs) << "\t"
           << percentage(qnameStats.NR, total_pairs) << "\t"
           << percentage(qnameStats.NN, total_pairs) << "\n";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Runner::update_pair_plots_from_records(const bam1_t* rec1, const bam1_t* rec2) {
    if (!rec1 || !rec2) return;
    if (rec1->core.tid != rec2->core.tid) return; // solo stesso riferimento

    const bam1_t* left = rec1;
    const bam1_t* right = rec2;

    if (left->core.pos > right->core.pos) {
        std::swap(left, right);
    }

    uint64_t dist = static_cast<uint64_t>(right->core.pos - left->core.pos);

    bool left_rev  = (left->core.flag  & BAM_FREVERSE);
    bool right_rev = (right->core.flag & BAM_FREVERSE);

    update_log_binned_distance(dist, graph.Ps_binned_dist_count);

    if      (!left_rev && right_rev)  update_log_binned_distance(dist, graph.fr_binned_dist_count);
    else if ( left_rev && !right_rev) update_log_binned_distance(dist, graph.rf_binned_dist_count);
    else if (!left_rev && !right_rev) update_log_binned_distance(dist, graph.ff_binned_dist_count);
    else if ( left_rev &&  right_rev) update_log_binned_distance(dist, graph.rr_binned_dist_count);
}

void Runner::update_reference_graph(const bam1_t* rec1, const bam1_t* rec2) {
    if (!rec1 || !rec2) return;

    uint32_t tid1 = rec1->core.tid;
    uint32_t tid2 = rec2->core.tid;

    if (tid1 == tid2) {
        ++graph.graph_intra_count[tid1];
    } else {
        uint32_t a = std::min(tid1, tid2);
        uint32_t b = std::max(tid1, tid2);
        ++graph.graph_inter_count[{a, b}];
    }
}

//se le read partono uguali ma finiscono diverse lo considero un walk in ogni caso, anche pairtools dovrebbe fare così, se è un miglioramento si può pensare ad un'implementazione
uint16_t Runner::Alignstarts(const bam1_t* b){//legge il cigar e riporta le basi segnate nel read a sinistra prima delle basi mappate
	const uint32_t* cigar = bam_get_cigar(b);
	std::size_t bases=0;

	for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        int   op = bam_cigar_op(cigar[i]);//cigar operator character
        int  len = bam_cigar_oplen(cigar[i]);//quante basi per lettera (es.50S)

		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {break;}//inizia l'allineamento con M,X,=
		if (op == BAM_CSOFT_CLIP || op == BAM_CINS) {bases+=len;}//sommo quante basi S e I ci sono prima che la read si allinei al riferimento
	}
	return bases;
}

int Runner::Alignend(const bam1_t* b) {
    const uint32_t* cigar = bam_get_cigar(b);
    int qlen = bam_cigar2qlen(b->core.n_cigar, cigar);

    int trailing = 0;
    for (int i = (int)b->core.n_cigar - 1; i >= 0; --i) {
        int op  = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) break;
        if (op == BAM_CSOFT_CLIP || op == BAM_CINS) trailing += len;
    }

    return qlen - trailing;
}

int Runner::inter_align_gap_on_query(const bam1_t* left_seg, const bam1_t* right_seg) {
    int left_end    = Alignend(left_seg);
    int right_start = Alignstarts(right_seg);
    return right_start - left_end;   // >0 gap, <0 overlap
}

void Runner::qname_stats(Bam_record_vector &group) {
    std::size_t begin = 0;
    std::size_t end = 0;

    //vettori di indici non di reads
    std::vector<std::size_t> r1_mapped; r1_mapped.reserve(5);
    std::vector<std::size_t> r2_mapped; r2_mapped.reserve(5); 

    std::vector<std::size_t> r1_all; r1_all.reserve(5);
    std::vector<std::size_t> r2_all; r2_all.reserve(5);

    std::size_t secondary_r1 = 0, secondary_r2 = 0;
    std::size_t primary_r1 = (std::size_t)-1, primary_r2 = (std::size_t)-1;
    std::size_t supplementary_r1 = 0, supplementary_r2 = 0;
    
	bool mapR1 = false;
	bool mapR2 = false;

    bool r1_unresolved = false;
    bool r2_unresolved = false;
    bool rescued = false;
    bool is_dupl = false;
    bool R1_chim = false;
    bool R2_chim = false;
    std::size_t inner = 0, outer = 0, other = 0;

	const std::size_t NO_INDEX = static_cast<std::size_t>(-1);
    std::size_t plot_r1 = NO_INDEX;
    std::size_t plot_r2 = NO_INDEX;

    std::size_t MAX_INTER_ALIGN_GAP=20;

    Maptype a = Maptype::N;
    Maptype b = Maptype::N;

    while (begin < group.size()) {

        const char* qname = bam_get_qname(group[begin]);
        end = begin + 1;
        while (end < group.size() && strcmp(bam_get_qname(group[end]), qname) == 0) ++end;

        ++pairStats.n_group;

        r1_mapped.clear();
        r2_mapped.clear();

        r1_all.clear();
        r2_all.clear();

        secondary_r1 = secondary_r2 = 0;
        supplementary_r1 = supplementary_r2 = 0;
        primary_r1 = primary_r2 = (std::size_t)-1;

		mapR1 = false;
		mapR2 = false;

        r1_unresolved = false;
        r2_unresolved = false;
        rescued = false;
        is_dupl = false;
        R1_chim = false;
        R2_chim = false;

    	plot_r1 = NO_INDEX;
    	plot_r2 = NO_INDEX;

        a = Maptype::N;
        b = Maptype::N;

        // -------------------------
        // 1) Raccolta dati
        // -------------------------
        //total_read = end - begin;

        for (std::size_t j = begin; j < end; ++j) {	

            auto flag = group[j]->core.flag;

	      

            if (flag & BAM_FSECONDARY) {

                if(flag & BAM_FREAD1) ++secondary_r1;
                if(flag & BAM_FREAD2) ++secondary_r2;

            } else { 

                if (flag & BAM_FREAD1) {
                    r1_all.push_back(j); //che sia mappata o non mappata

                    if (!(flag & BAM_FUNMAP)) {//se è mappata
                        r1_mapped.push_back(j);
                        if (flag & BAM_FSUPPLEMENTARY) ++supplementary_r1;
                        else {
				primary_r1 =j;
				if (flag & BAM_FDUP) is_dupl=true;
				 mapR1=true;  
			     }
                    }

                } else if (flag & BAM_FREAD2){
                    r2_all.push_back(j); //che sia mappata o non mappata

                    if (!(flag & BAM_FUNMAP)) {//se è mappata
                        r2_mapped.push_back(j);
                        if (flag & BAM_FSUPPLEMENTARY) ++supplementary_r2;
                        else {
				primary_r2 =j;
				if (flag & BAM_FDUP) is_dupl=true;
  				 mapR2=true;
			     }

                    }

                }
            }
        }

		if(is_dupl) {
			++qnameStats.DD;
			begin=end;
			continue;
		}


		if (mapR1 && mapR2) {
			if(is_dupl) ++pairStats.dupl;
			++pairStats.two_side_mapped;
			if (group[primary_r1]->core.tid = group[primary_r2]->core.tid) ++pairStats.cis; 
		}
		if (mapR1 != mapR2)	++pairStats.one_side;
		if (!mapR1 && !mapR2) ++pairStats.UNmapped;

		if(is_dupl) {
                	++qnameStats.DD;
                	begin=end;
        		continue;
        	}


        bool r1_has_null = (r1_all.size() > r1_mapped.size());
        bool r2_has_null = (r2_all.size() > r2_mapped.size());

        // -------------------------
        // 2) N/U/M
        // -------------------------

        // ----- lato R1 -----
        if (r1_mapped.size() == 0) {
            a = Maptype::N;
        }
	else if (secondary_r1 > 0) {
	    a = Maptype::M;
        }
	else if (r1_all.size() >= 2) {
            r1_unresolved = true;
        }
        else if (primary_r1 == (std::size_t)-1 ||
                group[primary_r1]->core.qual < 1) {
            a = Maptype::M;
        }
        else if (r1_mapped.size() == 1 && !r1_has_null && supplementary_r1 == 0) {
            a = Maptype::U;
        }
        else {
            // qui non decido ancora: potrebbe essere walk oppure rescue
            r1_unresolved = true;
        }

        // ----- lato R2 -----
        if (r2_mapped.size() == 0) {
            b = Maptype::N;
        }
	else if (secondary_r2 > 0) {
            b = Maptype::M;
        }
        else if (r2_all.size() >= 2) {
            r2_unresolved = true;
        }
        else if (primary_r2 == (std::size_t)-1 ||
                group[primary_r2]->core.qual < 1) {
            b = Maptype::M;
        }
        else if (r2_mapped.size() == 1 && !r2_has_null && supplementary_r2 == 0) {
            b = Maptype::U;
        }
        else {
            r2_unresolved = true;
        }

        // -------------------------
        // 3) WALK / RESCUE /CASI AMBIGUI
        // -------------------------

        if (r1_unresolved || r2_unresolved){
++qnameStats.dbg_unresolved;
           std::size_t total_all = r1_all.size() + r2_all.size();
		   const std::size_t NO_INDEX = static_cast<std::size_t>(-1);
			inner = outer = other = NO_INDEX;

            // se non è il classico 2+1, lo mando in WW
            if (total_all != 3) {
++qnameStats.dbg_not3;
                ++qnameStats.WW;
                begin = end;
                continue;
            }

            // pattern 2+1
            bool is_2plus1 =
                (r1_all.size() == 2 && r2_all.size() == 1) ||
                (r2_all.size() == 2 && r1_all.size() == 1);

            if (!is_2plus1) {
++qnameStats.dbg_not_2plus1;
                ++qnameStats.WW;
                begin = end;
                continue;

            } else {

				//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||			
				// -------------------------
				// identifica il lato chimerico
				// -------------------------
				bool chim_on_r1 = (r1_all.size() == 2); //così so cos'è e che uno è il contrario dell'altro
				R1_chim = chim_on_r1;
				R2_chim = !chim_on_r1;

				// riferimenti comodi
				auto& chim_all    = chim_on_r1 ? r1_all    : r2_all; //se lato chimerico R1 allora chim_all=R1_all se no R2_all
				auto& chim_mapped = chim_on_r1 ? r1_mapped : r2_mapped; //se lato chimerico ampped R1 allora chim_mapped=R1_mapped se no R2_mapped

				auto& other_all    = chim_on_r1 ? r2_all    : r1_all;
				auto& other_mapped = chim_on_r1 ? r2_mapped : r1_mapped;

				// -------------------------
				// scegli OUTER e INNER sul lato chimerico
				// -------------------------
				// Caso A: due mapped -> ho davvero inner e outer
				if (chim_mapped.size() == 2) {
					inner = chim_mapped[0];
					outer = chim_mapped[1];

					// ordino sulla query: inner è il più interno, outer il più esterno
					if (Alignstarts(group[inner]) > Alignstarts(group[outer])) {
						std::swap(inner, outer);
					}
				}
				// Caso B: un solo mapped -> quello è outer, inner è "assente/ignorabile"
				else if (chim_mapped.size() == 1) {
					outer = chim_mapped[0];
					inner = NO_INDEX;
				}
				// Caso C: nessun mapped sul lato chimerico -> non posso rescueare
				else {
++qnameStats.dbg_unmapped;
					++qnameStats.WW;
					begin = end;
					continue;
				}

				// -------------------------
				// scegli OTHER sul lato opposto
				// -------------------------
				// se il lato opposto non ha mapped, other resta NO_INDEX
				if (!other_mapped.empty()) {
					other = other_mapped[0];
				}
				//FINE INDENTIFICAZIONE LATO CHIMERICO
				//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

				Maptype other_type = R1_chim ? b : a;

				// 1) OUTER deve esistere ed essere buono
				bool outer_missing = (outer == NO_INDEX);
				bool outer_bad = outer_missing || (group[outer]->core.qual < 1);

				if (outer_bad) {
if(outer_missing){++qnameStats.dbg_outer_noindex;}
++qnameStats.dbg_outer_bad;
 if (other == NO_INDEX || other_type == Maptype::N) ++qnameStats.dbg_outer_bad_otherN;
    else if (other_type == Maptype::M) ++qnameStats.dbg_outer_bad_otherM;
    else if (other_type == Maptype::U) ++qnameStats.dbg_outer_bad_otherU;
					++qnameStats.WW;
					begin = end;
					continue;
				}

				// 2) controllo se INNER è ignorabile
				bool inner_missing = (inner == NO_INDEX);
				bool inner_bad = inner_missing || (group[inner]->core.qual < 1);

				// true se il gap/overlap rende inner non affidabile / ignorable
				bool inner_null_like = false;
				if (!inner_missing) {
					int qgap = inter_align_gap_on_query(group[inner], group[outer]);
					inner_null_like = (qgap <  MAX_INTER_ALIGN_GAP); //così prende il gap piccolo (e ache l'oevrlap)!!!!!!!!!!!!!!!!!!!!!!
 if(inner_null_like){++qnameStats.dbg_gap_large;}
				}

				bool ignore_inner = inner_bad || inner_null_like;

				// 3) se OTHER è N o M -> rescue diretto
				if (other == NO_INDEX || other_type == Maptype::N) {
++qnameStats.dbg_other_no_index;
++qnameStats.dbg_other_type_N;
					rescued = true;
					if (R1_chim) {
						a = Maptype::R;
						b = Maptype::N;
					} else {
						b = Maptype::R;
						a = Maptype::N;
					}
				}
				else if (other_type == Maptype::M) {
++qnameStats.dbg_other_type_M;
					rescued = true;
					if (R1_chim) {
						a = Maptype::R;
						b = Maptype::M;
					} else {
						b = Maptype::R;
						a = Maptype::M;
					}
				}

				// 4) se OTHER è U e INNER è ignorabile -> rescue diretto
				else if (other_type == Maptype::U && ignore_inner) {
++qnameStats.dbg_ignore_inner;
					rescued = true;
					if (R1_chim) a = Maptype::R;
					else         b = Maptype::R;
				}

				// 5) se OTHER è U e INNER è buono -> criteri geometrici
				else if (other_type == Maptype::U) {

					bool rev_i  = group[inner]->core.flag & BAM_FREVERSE;
					bool rev_ot = group[other]->core.flag & BAM_FREVERSE;

					bool cis =
						(group[inner]->core.tid == group[other]->core.tid);

					bool facing =
						((!rev_i && rev_ot && group[inner]->core.pos <= group[other]->core.pos) ||
						( rev_i && !rev_ot && group[other]->core.pos <= group[inner]->core.pos));

					bool distance =
						(llabs(group[inner]->core.pos - group[other]->core.pos) <= 2000);

					if (cis && facing && distance) {
++qnameStats.dbg_geom_pass;
						rescued = true;
						if (R1_chim) a = Maptype::R;
						else         b = Maptype::R;
					} else {
						++qnameStats.WW;
++qnameStats.dbg_geom_fail;
						begin = end;
						continue;
					}
				} else {// 6) fallback: qualunque caso non coperto -> WW
++qnameStats.fallback;
					++qnameStats.WW;
					begin = end;
					continue;
				}

			}

        }

		 if (a == Maptype::U && b == Maptype::U) {
            plot_r1 = primary_r1;
            plot_r2 = primary_r2;
        }
        else if (a == Maptype::R && b == Maptype::U) {
            plot_r1 = outer;   // lato R1 rescued
            plot_r2 = other;   // lato R2 unique
        }
        else if (a == Maptype::U && b == Maptype::R) {
            plot_r1 = other;   // lato R1 unique
            plot_r2 = outer;   // lato R2 rescued
        }

        if (plot_r1 != NO_INDEX && plot_r2 != NO_INDEX) {
            update_pair_plots_from_records(group[plot_r1], group[plot_r2]);
			update_reference_graph(group[plot_r1], group[plot_r2]);
        }
        // -------------------------
        // 6) classificazione finale
        // -------------------------
        bool classified = false;

		

        if (a == Maptype::U && b == Maptype::R) {
            ++qnameStats.UR;
            classified = true;
        }

        if (!classified) {

            if ((int)a < (int)b) std::swap(a, b);

            if (a == Maptype::U && b == Maptype::U) { ++qnameStats.UU; classified = true; }
            else if (a == Maptype::M && b == Maptype::M) { ++qnameStats.MM; classified = true; }
            else if (a == Maptype::N && b == Maptype::N) { ++qnameStats.NN; classified = true; }
            else if (a == Maptype::M && b == Maptype::U) { ++qnameStats.MU; classified = true; }
            else if (a == Maptype::U && b == Maptype::N) { ++qnameStats.NU; classified = true; }
            else if (a == Maptype::M && b == Maptype::N) { ++qnameStats.NM; classified = true; }
            else if (a == Maptype::R && b == Maptype::U) { ++qnameStats.RU; classified = true; }
            else if (a == Maptype::R && b == Maptype::M) { ++qnameStats.MR; classified = true; }
            else if (a == Maptype::R && b == Maptype::N) { ++qnameStats.NR; classified = true; }
        }

        if (!classified) {
            ++qnameStats.LOST;
        }

        begin = end;
    }
	}
}

double Runner::percentage(std::size_t value, double total) {
    if (total == 0.0) return 0.0;
    return 100.0 * static_cast<double>(value) / total;
}

long double Runner::update_mean_tlen(long double prev_mean,std::uint64_t k, bam1_t* bamdata){  //<x>
    long double xk = std::abs((long double)bamdata->core.isize);  // TEN dLel record
    return (xk / k) + ((k - 1) / (long double)k) * prev_mean;									
}

long double Runner::update_quadratic_mean_tlen(long double prev_qmean,std::uint64_t k, bam1_t* bamdata){ //<x^2> FORSE SBAGLIATA
	long double xk = std::abs((long double)bamdata->core.isize);  // TLEN del record
	long double xk2 = xk * xk;
 //   return (pow(xk,2) / k) + ((k - 1) / (long double)k) * pow(prev_mean,2);
	return prev_qmean + (xk2 - prev_qmean) / (long double)k;
}

double Runner::error_rate(uint64_t mismatched_bases,uint64_t total_base){
	return (total_base==0) ? 0.0 : (long double)mismatched_bases/total_base;
}

void Runner::histo_global_distance (std::unordered_map<uint64_t,uint64_t>& global_dist_count){
	std::fstream myfile;
	myfile.open("Pair_by_global_distance.txt",std::ios::out); //agginugere il path (i vecchi dati vengono cancellati e sovrascritti) 

	if(!myfile.is_open()){
		std::cout<<"pair_by_global_distance not open"<<std::endl;
		return;
	}

	myfile << "distance" << "\t" << "count" << "\n";
	for (auto i = global_dist_count.begin(); i != global_dist_count.end(); ++i) {
		myfile << i->first << "\t" << i->second<< "\n";
	}
	myfile.close();  
}

void Runner::histo_chrom_distance (std::map<uint32_t,std::unordered_map<uint64_t,uint64_t>>& chrom_dist_count) { 
	std::fstream myfile;
	myfile.open("Pair_chromosome_by_distance.txt",std::ios::out); //aggiungere il path (i vecchi dati vengono cancellati e sovrascritti)

	if(!myfile.is_open()){
		std::cout<<"pair_chromosome_by_distance not open"<<std::endl;
		return;
	}
	
	myfile << "\n# chromosome" <<"\t"<<"distance"<<"\t"<<"counter"<< "\n";


	for (auto i = chrom_dist_count.begin(); i != chrom_dist_count.end(); ++i) {

    	uint32_t chrom = i->first;
    	const auto& dist_map = i->second;

    	for (auto j = dist_map.begin(); j != dist_map.end(); ++j) {
       	 myfile <<i->first<< "\t" << j->first << "\t" << j->second << "\n";
    	}
	}
	myfile.close();  

}

void Runner::flag_inspector (bam1_t* bamdata) {
	uint16_t flag= bamdata-> core.flag;

	//sui singoli record

	if (flag & BAM_FQCFAIL){++readStats.qc_fail;} //return?
	if (flag & BAM_FUNMAP) {++readStats.unmapped;} //così i mapped sono di tutti come in samtools ma ho la richiesta una volta sola e non dentro e fuori dall'else
	if (flag & BAM_FPROPER_PAIR) {++pairStats.proper_pairs;}

	if(flag & BAM_FSUPPLEMENTARY) {
		++readStats.supplementary;
		return;  // o si tolgono i return e mettere le cose di coppia nell'else come fa samtools flagstats
	} else if (flag & BAM_FSECONDARY) {
		++readStats.secondary;
		return;
	} else {readStats.primary++;}

	//sulle coppie, faccio tutte e due insieme così ho già filtrato dall calcolo le coppie supplementary e secondary
	
	if(flag & BAM_FPAIRED || flag & BAM_FPROPER_PAIR) {
		++pairStats.pairN;

		if(flag & BAM_FDUP) {++pairStats.duplicated;} //WARNING! se c'è solo una read1 ma segnata come duplicato così non la conto nelle rad (basta toglier l'else if)
		else if (flag & BAM_FREAD1) { // così ne prendo solo una e non due non so se ha senso
			++pairStats.read1;

			//if ((flag & BAM_FUNMAP && !(flag & BAM_FMUNMAP))^(flag & BAM_FMUNMAP && !(flag & BAM_FUNMAP))) {++pairStats.UMone_sided;} // statistica fatta sul singolo se no è doppia
			//else if (flag & BAM_FUNMAP && (flag & BAM_FMUNMAP)) {++pairStats.UNmapped;}
	
			if(!(flag & BAM_FUNMAP) && !(flag & BAM_FMUNMAP)) {
				//++pairStats.UMtwo_sided;
				pairStats.good_read1=true;
				++pairStats.good_pairs;
			}
		}
		else if (flag & BAM_FREAD2) {
			++pairStats.read2;
			if(!(flag & BAM_FUNMAP) && !(flag & BAM_FMUNMAP)) {
				pairStats.good_read2=true;
			}
		}	
	} 
}
	

void Runner::processReads(Bam_record_vector &vectorbox) {
	if(userInput.single_read_stats || userInput.hist_global || userInput.hist_by_chrom){

		uint32_t chrom=0;
		uint64_t dist=0;

		if(userInput.single_read_stats){ 
			for(int i=0;i<vectorbox.size();++i){
				pairStats.good_read1=false;
				pairStats.good_read2=false;
				++readStats.readN;
	
				if (!(vectorbox[i]->core.flag & BAM_FUNMAP) && !(vectorbox[i]->core.flag & BAM_FSUPPLEMENTARY) && !(vectorbox[i]->core.flag & BAM_FSECONDARY) && vectorbox[i]->core.qual==0) {++readStats.mapQ0;}
				flag_inspector(vectorbox[i]);
 
				if (pairStats.good_read1 && vectorbox[i]->core.tid == vectorbox[i]->core.mtid) {
					//++pairStats.sameCr;
					if(((vectorbox[i]->core.flag & BAM_FREVERSE) != (vectorbox[i]->core.flag & BAM_FMREVERSE)) && std::abs((long double)vectorbox[i]->core.isize)>0){//così hanno sempre orientamenti opposti
						if(std::abs((long double)vectorbox[i]->core.isize)<10000){
							++readStats.avf_counter;
							readStats.mean_insert_filtr = update_mean_tlen(readStats.mean_insert_filtr, readStats.avf_counter, vectorbox[i]);
							readStats.quadratic_mean_filtr=update_quadratic_mean_tlen(readStats.quadratic_mean_filtr,readStats.avf_counter, vectorbox[i]);
							}
						++readStats.av_counter;			
						readStats.mean_insert = update_mean_tlen(readStats.mean_insert, readStats.av_counter, vectorbox[i]);   
						//	readStats.quadratic_mean=update_quadratic_mean_tlen(readStats.mean_insert,av_counter, bamdata);
						readStats.quadratic_mean=update_quadratic_mean_tlen(readStats.quadratic_mean,readStats.av_counter, vectorbox[i]);
						
					}

					dist=llabs(vectorbox[i]->core.pos - vectorbox[i]->core.mpos); // dovrebbero essere degli uint64_t quindi non serve forzare il double ne arrotondare
					
					if(userInput.hist_global){ //HISTO_GLOBAL_DATA
						++global_dist_count[dist]; 
					}
					
					if(userInput.hist_by_chrom){	 //HISTO_CHROM_DATA
						chrom=vectorbox[i]->core.tid;
						++chrom_dist_count[chrom][dist];
					}
						
					/*if(pairStats.good_read1 ){ 
						bool r1_rev = (vectorbox[i]->core.flag & BAM_FREVERSE);//POTREI USARLI IN TUTTA
						bool r2_rev = (vectorbox[i]->core.flag & BAM_FMREVERSE);

						update_log_binned_distance(dist,graph.Ps_binned_dist_count); 

						if      (!r1_rev && r2_rev)  update_log_binned_distance(dist, graph.fr_binned_dist_count);// Forward-Reverse (Inward) ---> <---
						else if (r1_rev && !r2_rev)  update_log_binned_distance(dist, graph.rf_binned_dist_count);// Reverse-Forward (Outward) <--- --->
						else if (!r1_rev && !r2_rev) update_log_binned_distance(dist, graph.ff_binned_dist_count);// Forward-Forward (Same-strand) ----> --->
						else if (r1_rev && r2_rev) 	 update_log_binned_distance(dist, graph.rr_binned_dist_count);// Reverse-Reverse (Same-strand) <---- <----
											
					}*/

					if (pairStats.good_read1 || pairStats.good_read2) { 
			      	// if (!(vectorbox[i]->core.flag & BAM_FUNMAP) &&  !(vectorbox[i]->core.flag & BAM_FSECONDARY) && !(vectorbox[i]->core.flag & BAM_FSUPPLEMENTARY)){
                        uint8_t* nm_ptr = bam_aux_get(vectorbox[i], "NM");//diff tra a read e il riferimento
						uint64_t nm = nm_ptr ? bam_aux2i(nm_ptr) : 0;

						readStats.mismatched_bases += nm;  
						uint64_t aligned = bam_cigar2rlen(vectorbox[i]->core.n_cigar, bam_get_cigar(vectorbox[i])); //bam_cigar2rlen(int n_cigar, const uint32_t *cigar):This function returns the sum of the lengths of the M, I, S, = and X operations in @p cigar (these are the operations that "consume" query bases
						readStats.total_mapped_base += aligned;
					} 
				}

				if ((readStats.readN % 2) == 0) {
					pairStats.pairN=readStats.readN/2;
				} else {pairStats.pairN += (readStats.readN % 2) * 2 >= 2 ? 1 : 0;} //arromtonda all'intero più vicino 
				//readStats.error_rate=error_rate(mismatched_bases,total_base);
			}
		}

		if(userInput.pair_read_stats) qname_stats(vectorbox);
	}
}


void Runner::output(){
		std::cout<<"SINGLE READ STATISTICS:"<<std::endl;
		std::cout<<"Tot_raw_record: "<<readStats.readN<<std::endl;
		std::cout<<"Non_primary: "<<readStats.secondary+readStats.supplementary<<std::endl;
		std::cout<<"Primary_read: "<<readStats.readN-(readStats.secondary+readStats.supplementary)<<std::endl;
		std::cout<<"Supplementary: "<<readStats.supplementary<<std::endl; 
		std::cout<<"Read1: "<<pairStats.read1<<std::endl;
		std::cout<<"Read2: "<<pairStats.read2<<std::endl;
		std::cout<<"Reads_mapped: "<<readStats.readN-readStats.unmapped<<std::endl;
		std::cout<<"Unmapped: "<<readStats.unmapped<<std::endl;        
		std::cout<<"Proper_pairs: "<<pairStats.proper_pairs<<std::endl;;
		std::cout<<"Record_duplicated: "<<pairStats.duplicated<<std::endl;
		std::cout<<"MapQ0: "<<readStats.mapQ0<<std::endl;
		std::cout<<"Qc_fail: "<<readStats.qc_fail<<std::endl; 
		std::cout<<"Pairs: "<<pairStats.pairN<<std::endl;
		std::cout<<"insert_size_avarage_filtred: "<<readStats.mean_insert_filtr<<std::endl;
		std::cout<<"SD_filtr: "<<sqrt(readStats.quadratic_mean_filtr-pow(readStats.mean_insert_filtr,2))<<std::endl; //rad(<x^2>-<x>^2)
		std ::cout<<"insert_size_average: "<<readStats.mean_insert<<std::endl;
		std::cout<<"SD: "<<std::sqrt(readStats.quadratic_mean -(readStats.mean_insert * readStats.mean_insert))<<std::endl;
		//std ::cout<<"SD_insert_size:"<<readStats.sd_insert<<std::endl;
		//std::cout << "insert size average:\t" << readStats.mean_insert << std::endl;
		std::cout<<"error_rate: "<<error_rate(readStats.mismatched_bases,readStats.total_mapped_base)<<std::endl;
		std::cout<<"|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<std::endl;
		std::cout<<"PAIR READS STATISTICS:"<<std::endl;
		std::cout<<"Pairs_two_sided_mapped: "<<pairStats.two_side_mapped<<std::endl;
		std::cout<<"%Two_sided_mapped: "<< ((pairStats.two_side_mapped*100)/(long double)pairStats.pairN)<<std::endl;
                std::cout<<"Pairs_one_sided_mapped: "<<pairStats.one_side<<std::endl;
		std::cout<<"%One_sided: "<< ((pairStats.one_side*100)/(long double)pairStats.pairN)<<std::endl;
                std::cout<<"Pairs_unmapped: "<<pairStats.UNmapped<<std::endl;
		std::cout<<"%Unmapped: "<<((pairStats.UNmapped*100)/(long double)pairStats.pairN)<<std::endl;
                std::cout<<"CIS: "<<pairStats.cis<<std::endl;
		std::cout<<"Primary_dupl "<<pairStats.dupl<<std::endl;
		std::cout<<"unique_primary: "<<pairStats.two_side_mapped-pairStats.dupl<<std::endl;//mapped primary after deduplication
		std::cout<<std::endl;
		std::cout<<"%Two_sided_mapped: "<< ((pairStats.two_side_mapped*100)/(long double)pairStats.pairN)<<std::endl;
		std::cout<<"%One_sided: "<< ((pairStats.one_side*100)/(long double)pairStats.pairN)<<std::endl;
		std::cout<<"%Unmapped: "<<((pairStats.UNmapped*100)/(long double)pairStats.pairN)<<std::endl;
		std::cout<<"%CIS: "<<((pairStats.cis*100)/(long double)pairStats.two_side_mapped)<<std::endl;
		std::cout<<"%primary_dupl: "<<((pairStats.dupl*100)/(long double)pairStats.two_side_mapped)<<std::endl;
		std::cout<<std::endl;
		std::cout<<"UU"<<":"<<"MM"<<":"<<"NN"<<":"<<"UM"<<":"<<"UN"<<":"<<"NM"<<"\t"<<qnameStats.UU<<":"<<qnameStats.MM<<":"<<qnameStats.NN<<":"<<qnameStats.MU<<":"<<qnameStats.NU<<":"<<qnameStats.NM<<std::endl;
		std::cout<<"DD"<<":"<<"WW"<<":"<<"UR"<<":"<<"RU"<<":"<<"RN"<<":"<<"RM"<<"\t"<<qnameStats.DD<<":"<<qnameStats.WW<<":"<<qnameStats.UR<<":"<<qnameStats.RU<<":"<<qnameStats.NR<<":"<<qnameStats.MR<<std::endl;

		std::size_t  tot_qname_stats= qnameStats.UU+qnameStats.MM+qnameStats.NN+qnameStats.MU+qnameStats.NU+qnameStats.NM+qnameStats.DD+qnameStats.WW+qnameStats.UR+qnameStats.RU+qnameStats.NR+qnameStats.MR;
		std::cout<<"%two_sided_mapped"<< (qnameStats.UU+qnameStats.UR+qnameStats.RU)/(long double)tot_qname_stats<<std::endl;
		std::cout<<"two_sided_mapped"<< qnameStats.UU+qnameStats.UR+qnameStats.RU<<std::endl;
		std::cout<<"total read"<<tot_qname_stats<<std::endl;
		//std::cout<<"funzione di overlapping"<<pairStats.counter1<<std::endl;

		std::cout << "unresolved: " << qnameStats.dbg_unresolved << "\n";
		std::cout << "not_2plus1->WW: " << qnameStats.dbg_not_2plus1 << "\n";
		std::cout << "not 3"<<  qnameStats.dbg_not3 << std::endl;
		std::cout << "no mapped chim side"<<  qnameStats.dbg_unmapped << std::endl;
		std::cout << "outer_bad->WW: " << qnameStats.dbg_outer_bad << "\n";
		std::cout << "outer noindex (not possible)"<<  qnameStats.dbg_outer_noindex << std::endl;
		std::cout << "outer bad other M"<<  qnameStats.dbg_outer_bad_otherM << std::endl;
		std::cout << "outer bad other U"<<  qnameStats.dbg_outer_bad_otherU << std::endl;
		std::cout << "outer bad other N"<<  qnameStats.dbg_outer_bad_otherN << std::endl;
		std::cout << "other==NO_INDEX: " << qnameStats.dbg_other_no_index << "\n";
		std::cout << "other_type==N: " << qnameStats.dbg_other_type_N << "\n";
		std::cout << "other_type==M: " << qnameStats.dbg_other_type_M << "\n";
		std::cout << "ignore_inner->R: " << qnameStats.dbg_ignore_inner << "\n";
		std::cout << "ignore inner per gap<20"<<  qnameStats.dbg_gap_large << std::endl;
		std::cout << "geom pass->R: " << qnameStats.dbg_geom_pass << "\n";
		std::cout << "geom fail->WW: " << qnameStats.dbg_geom_fail << "\n";
		std::cout << "fall back: " << qnameStats.fallback << "\n";

std::cout<<"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"<<std::endl;
		std::cout<<"mapped:"<<(readStats.readN-(readStats.secondary+readStats.supplementary))-readStats.unmapped<<std::endl;
                std::cout<<"%mapped:"<< 100-((readStats.unmapped*100)/(long double)readStats.readN-(readStats.secondary+readStats.supplementary))<<"%"<<std::endl; 
                std::cout<<"%MapQ0:"<< ((readStats.mapQ0*100)/(long double)(readStats.readN-(readStats.secondary+readStats.supplementary)))<<"%"<<std::endl;
                std::cout<<"MapQ0:"<<readStats.mapQ0<<std::endl;
                std::cout<<"Pairs:"<<(readStats.readN-(readStats.secondary+readStats.supplementary))/2<<std::endl;
                std::cout<<"Pairs_tot_veri?"<<(readStats.readN)/2<<std::endl; 

	}

void Runner::data_vector(Bam_record_vector &vectorbox,samFile *fp_in,bam_hdr_t *bamHdr){
	vectorbox.clear();
	for (int i=0; i<vectorbox.get_size_wanted();++i){ 
		if(!vectorbox.add_record(fp_in,bamHdr)) break;
	}
}

void Runner::data_vector(Bam_record_vector &vectorbox, bam1_t *bridge_read,bool &first, samFile *fp_in, bam_hdr_t *bamHdr){
	const char* qname;
	const char* current_qname;
	bool bridge=true;
	vectorbox.clear();

	for (int i=0;i<vectorbox.get_size_wanted();++i){
		if(bridge){
			if(first){
				vectorbox.add_record(fp_in,bamHdr);
				first =false;
			}else{
				vectorbox.push_back(bridge_read);
			}
			bridge=false;
			continue;
		}
		vectorbox.add_record(fp_in,bamHdr);
	}

	qname=bam_get_qname(vectorbox[vectorbox.size()-1]);
	for(;;){
		if(sam_read1(fp_in, bamHdr, bridge_read)>=0){ 
			current_qname=bam_get_qname(bridge_read);
			if(strcmp(current_qname, qname) != 0){
				break;}
			vectorbox.push_back(bridge_read); 
		}else{break;}                 
	}
}



void Runner::run() {

	std::size_t numFiles = userInput.inFiles.size();
	lg.verbose("Processing " + std::to_string(numFiles) + " files");
	
	for (uint32_t i = 0; i < numFiles; ++i) {

		global_dist_count.clear();//svuotare le mappe prima di ogni file o si mescoleranno (se non è l'obiettivo)s
		chrom_dist_count.clear();
		graph.Ps_binned_dist_count.clear();//!!!!!!!!!!!non volgio sia comulativo tra file diversi COME FACCIO SE OGNI FILE OUTPUT MI SOVRASCRIVE QUELLO PERIMA PER PIÙ FILE
		graph.ff_binned_dist_count.clear();
		graph.fr_binned_dist_count.clear();
		graph.rf_binned_dist_count.clear();
		graph.rr_binned_dist_count.clear();

		std::string file = userInput.file('r', i);
		std::string ext = getFileExt(file);
		
		samFile *fp_in = hts_open(userInput.file('r', i).c_str(),"r"); 
		if (!fp_in) {std::cout<<"hts_open has failed"<<std::endl;}
		bam_hdr_t *bamHdr = sam_hdr_read(fp_in); 
		if (!bamHdr) {std::cout<<"sam_hdr_read has failed"<<std::endl;}
		
		htsThreadPool tpool_read = {NULL, 0};
		tpool_read.pool = hts_tpool_init(userInput.decompression_threads);
		if (tpool_read.pool) {	hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &tpool_read);
		} else { lg.verbose("Failed to generate decompression threadpool with " + std::to_string(userInput.decompression_threads) + " threads. Continuing single-threaded");}

		bool qname_sorted =(std::string(bamHdr->text, bamHdr->l_text).find("SO:queryname") != std::string::npos); // perchè la funzione string.find() ritorna npos;
		if (!userInput.single_read_stats && !userInput.pair_read_stats) {//default
			userInput.single_read_stats = true;
    		if (qname_sorted) {
        	userInput.pair_read_stats = true;
    		} else { std::cout<<"Warning: input BAM file is not qname sorted, pair read statistics will not be computed."<<std::endl;}
		}
		if(userInput.pair_read_stats && !qname_sorted){
			userInput.pair_read_stats=false;
			std::cout<<"Error: to compute pair read statistics the input BAM file must be qname sorted."<<std::endl;
			exit(1);
		}
		std::size_t j=10000;//set real capacity
		bool first=true;
 
		Bam_record_vector records_vector(j); 
		bam1_t *bridge_read=bam_init1(); //È UN PUNTATORE

		while(!(records_vector.is_file_end())){ 
			if(userInput.pair_read_stats){
				data_vector(records_vector, bridge_read,first, fp_in, bamHdr);
				if((strcmp(bam_get_qname(records_vector[0]), bam_get_qname(records_vector[1])) != 0)) std::cout<<"!"<<std::endl;
			}else if(userInput.single_read_stats){
				data_vector(records_vector,fp_in,bamHdr);
			}
			processReads(records_vector);
    	}

//		estimate_insert_stats();
		if(userInput.hist_global){histo_global_distance(global_dist_count);}
		if(userInput.hist_by_chrom){histo_chrom_distance(chrom_dist_count);}

    		
		write_all_stats_file("all_stats.tsv",bamHdr);//!!!!!!!!!!!!!!!!!!!!!!!!!!!
		output();
		
		bam_hdr_destroy(bamHdr);
		bam_destroy1(bridge_read);
		sam_close(fp_in);
		if (tpool_read.pool)
		hts_tpool_destroy(tpool_read.pool);
	}
}


//////////////////////////////////////////////////////////////////////////////////////////class functions definition

Bam_record_vector::Bam_record_vector(std::size_t initial_capacity)
    : used(0), hiwater_data(0), size_wanted(0), file_end(false)
{
    slots.reserve(initial_capacity);
	size_wanted=initial_capacity;
    for (std::size_t i = 0; i < initial_capacity; ++i){ 
        slots.push_back(bam_init1());
	}
}

Bam_record_vector::~Bam_record_vector() {
   for (auto* b : slots) bam_destroy1(b);
}

Bam_record_vector::Bam_record_vector(Bam_record_vector&& other) noexcept //move contructor
    : slots(std::move(other.slots)),
    used(other.used),
    hiwater_data(other.hiwater_data),
	size_wanted(other.size_wanted),
    file_end(other.file_end)
{
	other.slots.clear();
    other.used= 0;
    other.hiwater_data= 0;
	other.size_wanted = 0;
    other.file_end = false;
}
Bam_record_vector& Bam_record_vector::operator=(Bam_record_vector&& other) noexcept{ //omve assignment operator
    if (this == &other) return *this;
    	for (auto* b : slots) bam_destroy1(b);

    	slots= std::move(other.slots);
    	used= other.used;
   		hiwater_data= other.hiwater_data;
		size_wanted = other.size_wanted;
    	file_end = other.file_end;

		other.slots.clear();
    	other.used= 0;
    	other.hiwater_data= 0;
		other.size_wanted = 0;
        other.file_end = false;
    return *this;
}

void Bam_record_vector::clear() noexcept { used=0;}
std::size_t Bam_record_vector::size() const noexcept { return used; }
std::size_t Bam_record_vector::capacity() const noexcept { return slots.size(); }
std::size_t Bam_record_vector::get_size_wanted() const noexcept {return size_wanted;}
bool Bam_record_vector::is_file_end() const noexcept {return file_end;}

bam1_t* Bam_record_vector::operator[](std::size_t i) noexcept { return slots.at(i); }
const bam1_t* Bam_record_vector::operator[](std::size_t i) const noexcept { return slots[i]; }

bool Bam_record_vector::add_record(samFile *fp_in,bam_hdr_t *bamHdr){
	if(used==slots.size()){expand(slots.empty() ? 10 : slots.size() * 2);}
	if(sam_read1(fp_in, bamHdr, slots[used])>=0){
		++used;
		return true;
	} else {
		file_end=true;
		return false;
	}
}

bam1_t* Bam_record_vector::push_back(const bam1_t* src) { // da sorgente al primo slot libero del vectorbox.
    if (used== slots.size())
        expand(slots.empty() ? 10 : slots.size() * 2);//se p empty per qualche motivo a cosa lo metto aiuto ahahha

    bam1_t* dst = slots[used];

	if (src->l_data > dst->m_data) {
    	bam1_t* new_dst = bam_dup1(src);
	 	if (!new_dst)
        throw std::bad_alloc();

        bam_destroy1(dst);
        slots[used] = new_dst;
        dst = new_dst;
    } else {
        bam_copy1(dst, src);
    }
    ++used;

    hiwater_data= std::max<int>(hiwater_data, dst->l_data);

	return dst;
}

void Bam_record_vector::expand(std::size_t new_capacity) {
    if (new_capacity <= slots.size()) return;
    slots.reserve(new_capacity);
    while (slots.size() < new_capacity) {
        bam1_t* b = bam_init1();
        if (!b) throw std::bad_alloc();

		if (hiwater_data > 0) {
   			b->data = (uint8_t*)malloc(hiwater_data);
   	 		if (!b->data) {
       			bam_destroy1(b);
       			throw std::bad_alloc();
   			}
   			b->m_data = hiwater_data;
    		b->l_data = 0;
		}
		slots.push_back(b);
    }
}

/*bool is_overlapping(bam1_t *a, bam1_t *b, int32_t max_overlap = 20) { 
    auto get_query_bounds = [](bam1_t *rec, int32_t &q_start, int32_t &q_end) {
        uint32_t *cigar = bam_get_cigar(rec);
        q_start = 0;
        q_end = 0;
        
        // Il "Soft Clipping" (op 4) all'inizio ci dice dove inizia il match sulla read
        if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
            q_start = bam_cigar_oplen(cigar[0]);
        }
        
        // Calcoliamo la fine sommando la lunghezza consumata sulla query dal CIGAR
        int32_t query_consumed = 0;
        for (uint32_t i = 0; i < rec->core.n_cigar; ++i) {
            if (bam_cigar_type(bam_cigar_op(cigar[i])) & 1) { // Consuma query (M, I, S)
                query_consumed += bam_cigar_oplen(cigar[i]);
            }
        }
        q_end = q_start + query_consumed;
    };

    int32_t start_a, end_a, start_b, end_b;
    get_query_bounds(a, start_a, end_a);
    get_query_bounds(b, start_b, end_b);

    // Calcolo dell'intersezione tra i due intervalli [start, end]
    int32_t overlap_start = std::max(start_a, start_b);
    int32_t overlap_end = std::min(end_a, end_b);
    int32_t overlap_len = overlap_end - overlap_start;

    // Se l'overlap è maggiore della soglia, sono sovrapposti (quindi NON è una chimera pulita)
    return (overlap_len > max_overlap);
}*/
