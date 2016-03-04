
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Sergey Koren beginning on 2016-FEB-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

/*
 * =====================================================================================
 *
 *       Filename:  kmer_count.c
 *
 *    Description:
 *
 *        Version:  0.1
 *        Created:  07/20/2013 17:00:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jason Chin,
 *        Company:
 *
 * =====================================================================================

 #################################################################################$$
 # Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
 #
 # All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted (subject to the limitations in the
 # disclaimer below) provided that the following conditions are met:
 #
 #  * Redistributions of source code must retain the above copyright
 #  notice, this list of conditions and the following disclaimer.
 #
 #  * Redistributions in binary form must reproduce the above
 #  copyright notice, this list of conditions and the following
 #  disclaimer in the documentation and/or other materials provided
 #  with the distribution.
 #
 #  * Neither the name of Pacific Biosciences nor the names of its
 #  contributors may be used to endorse or promote products derived
 #  from this software without specific prior written permission.
 #
 # NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
 # GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
 # BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 # WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 # USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 # OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 # SUCH DAMAGE.
 #################################################################################$$
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "falcon.H"

namespace FConsensus {
const uint32 KMERMATCHINC = 10000;

int compare_seq_coor(const void * a, const void * b) {
    const seq_coor_t * arg1 = (const seq_coor_t *)a;
    const seq_coor_t * arg2 = (const seq_coor_t *)b;
    return  (* arg1) - (* arg2);
}


kmer_lookup * allocate_kmer_lookup ( seq_coor_t size ) {
    kmer_lookup * kl;

    kl = (kmer_lookup *)  malloc( size * sizeof(kmer_lookup) );
    init_kmer_lookup( kl, size);
    return kl;
}

void init_kmer_lookup ( kmer_lookup * kl,  seq_coor_t size ) {
    seq_coor_t i;
    for (i=0; i<size; i++) {
        kl[i].start = INT_MAX;
        kl[i].last = INT_MAX;
        kl[i].count = 0;
    }
}


void free_kmer_lookup( kmer_lookup *  ptr) {
    free(ptr);
}

seq_array allocate_seq(seq_coor_t size) {
    seq_array sa;
    sa  = (seq_array) malloc( size * sizeof(base) );
    init_seq_array( sa, size);
    return sa;
}

void init_seq_array( seq_array sa, seq_coor_t size) {
    seq_coor_t i;
    for (i=0; i<size; i++) {
        sa[i] = 0xff;
    }
}

void free_seq_array( seq_array sa) {
    free(sa);
}

seq_addr_array allocate_seq_addr(seq_coor_t size) {
    return (seq_addr_array) calloc( size, sizeof(seq_addr));
}

void free_seq_addr_array(seq_addr_array sda) {
    free(sda);
}

seq_coor_t get_kmer_bitvector(seq_array sa, uint32 K) {
    uint32 i;
    seq_coor_t kmer_bv = 0;
    seq_coor_t kmer_mask;

    kmer_mask = 0;
    for (i = 0; i < K; i++) {
        kmer_mask <<= 2;
        kmer_mask |= 0x00000003;
    }

    for (i = 0; i < K; i++) {
        kmer_bv <<= 2;
        kmer_bv |= (uint32) sa[i];
    }

    return kmer_bv;
}

void add_sequence ( seq_coor_t start,
                    uint32 K,
                    const char * seq,
                    seq_coor_t seq_len,
                    seq_addr_array sda,
                    seq_array sa,
                    kmer_lookup * lk ) {

    seq_coor_t i;
    seq_coor_t kmer_bv;
    seq_coor_t kmer_mask;

    kmer_mask = 0;
    for (i = 0; i < K; i++) {
        kmer_mask <<= 2;
        kmer_mask |= 0x00000003;
    }

    for (i = 0; i < seq_len; i++) {
        switch ( seq[i] ) {
            case 'A':
                sa[ start + i ] = 0;
                break;
            case 'C':
                sa[ start + i ] = 1;
                break;
            case 'G':
                sa[ start + i ] = 2;
                break;
            case 'T':
                sa[ start + i ] = 3;
                break;
            default:
                sa[ start + i ] = 0;
        }
    }
    kmer_bv = get_kmer_bitvector( sa + start, K);
    for (i = 0; i < seq_len - K;  i++) {
        //fprintf(stderr, "%lu %lu\n", i, kmer_bv);
        //fprintf(stderr, "lk before init: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
        if (lk[kmer_bv].start == INT_MAX) {
            lk[kmer_bv].start = start + i;
            lk[kmer_bv].last = start + i;
            lk[kmer_bv].count += 1;
            //fprintf(stderr, "lk init: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
        } else {
            sda[ lk[kmer_bv].last ] = start + i;
            lk[kmer_bv].count += 1;
            lk[kmer_bv].last = start + i;
            //fprintf(stderr, "lk change: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
        }
        kmer_bv <<= 2;
        kmer_bv |= sa[ start + i + K];
        kmer_bv &= kmer_mask;
    }
}


void mask_k_mer(seq_coor_t size, kmer_lookup * kl, seq_coor_t threshold) {
    seq_coor_t i;
    for (i=0; i<size; i++) {
        if (kl[i].count > threshold) {
            kl[i].start = INT_MAX;
            kl[i].last = INT_MAX;
            //kl[i].count = 0;
        }
    }
}


kmer_match * find_kmer_pos_for_seq( const char * seq, seq_coor_t seq_len, uint32 K,
                    seq_addr_array sda,
                    kmer_lookup * lk) {
    seq_coor_t i;
    seq_coor_t kmer_bv;
    seq_coor_t kmer_mask;
    seq_coor_t kmer_pos;
    seq_coor_t next_kmer_pos;
    uint32 half_K;
    seq_coor_t kmer_match_rtn_allocation_size = KMERMATCHINC;
    kmer_match * kmer_match_rtn;
    base * sa;

    kmer_match_rtn = (kmer_match *) malloc( sizeof(kmer_match) );
    kmer_match_rtn->count = 0;
    kmer_match_rtn->query_pos = (seq_coor_t *) calloc( kmer_match_rtn_allocation_size, sizeof( seq_coor_t ) );
    kmer_match_rtn->target_pos = (seq_coor_t *) calloc( kmer_match_rtn_allocation_size, sizeof( seq_coor_t ) );

    sa = (base *)calloc( seq_len, sizeof(base) );

    kmer_mask = 0;
    for (i = 0; i < K; i++) {
        kmer_mask <<= 2;
        kmer_mask |= 0x00000003;
    }

    for (i = 0; i < seq_len; i++) {
        switch ( seq[i] ) {
            case 'A':
                sa[ i ] = 0;
                break;
            case 'C':
                sa[ i ] = 1;
                break;
            case 'G':
                sa[ i ] = 2;
                break;
            case 'T':
                sa[ i ] = 3;
                break;
             default:
                sa[ i ] = 0;
        }
    }


    kmer_bv = get_kmer_bitvector(sa, K);
    half_K = K >> 1;
    for (i = 0; i < seq_len - K;  i += half_K) {
        kmer_bv = get_kmer_bitvector(sa + i, K);
        if (lk[kmer_bv].start == INT_MAX) {  //for high count k-mers
            continue;
        }
        kmer_pos = lk[ kmer_bv ].start;
        next_kmer_pos = sda[ kmer_pos ];
        kmer_match_rtn->query_pos[ kmer_match_rtn->count ] = i;
        kmer_match_rtn->target_pos[ kmer_match_rtn->count ] = kmer_pos;
        kmer_match_rtn->count += 1;
        if (kmer_match_rtn->count > kmer_match_rtn_allocation_size - 1000) {
            kmer_match_rtn_allocation_size += KMERMATCHINC;
            kmer_match_rtn->query_pos = (seq_coor_t *) realloc( kmer_match_rtn->query_pos,
                                                                   kmer_match_rtn_allocation_size  * sizeof(seq_coor_t) );
            kmer_match_rtn->target_pos = (seq_coor_t *) realloc( kmer_match_rtn->target_pos,
                                                                    kmer_match_rtn_allocation_size  * sizeof(seq_coor_t) );
        }
        while ( next_kmer_pos > kmer_pos ){
            kmer_pos = next_kmer_pos;
            next_kmer_pos = sda[ kmer_pos ];
            kmer_match_rtn->query_pos[ kmer_match_rtn->count ] = i;
            kmer_match_rtn->target_pos[ kmer_match_rtn->count ] = kmer_pos;
            kmer_match_rtn->count += 1;
            if (kmer_match_rtn->count > kmer_match_rtn_allocation_size - 1000) {
                kmer_match_rtn_allocation_size += KMERMATCHINC;
                kmer_match_rtn->query_pos = (seq_coor_t *) realloc( kmer_match_rtn->query_pos,
                                                                       kmer_match_rtn_allocation_size  * sizeof(seq_coor_t) );
                kmer_match_rtn->target_pos = (seq_coor_t *) realloc( kmer_match_rtn->target_pos,
                                                                        kmer_match_rtn_allocation_size  * sizeof(seq_coor_t) );
            }
        }
    }
    free(sa);
    return kmer_match_rtn;
}

void free_kmer_match( kmer_match * ptr) {
    free(ptr->query_pos);
    free(ptr->target_pos);
    free(ptr);
}

aln_range* find_best_aln_range(kmer_match * km_ptr,
                              seq_coor_t K,
                              seq_coor_t bin_size,
                              seq_coor_t count_th) {
    seq_coor_t i;
    seq_coor_t j;
    seq_coor_t q_min, q_max, t_min, t_max;
    seq_coor_t * d_count;
    seq_coor_t * q_coor;
    seq_coor_t * t_coor;
    aln_range * arange;

    long int d, d_min, d_max;
    long int cur_score;
    long int max_score;
    long int max_k_mer_count;
    long int max_k_mer_bin;
    seq_coor_t cur_start;

    arange = (aln_range *)calloc(1 , sizeof(aln_range));

    q_min = INT_MAX;
    q_max = 0;
    t_min = INT_MAX;
    t_max = 0;

    d_min = INT_MAX;
    d_max = LONG_MIN;

    for (i = 0; i <  km_ptr->count; i++ ) {
        if ( km_ptr -> query_pos[i] < q_min) {
            q_min =  km_ptr->query_pos[i];
        }
        if ( km_ptr -> query_pos[i] > q_max) {
            q_max =  km_ptr->query_pos[i];
        }
        if ( km_ptr -> target_pos[i] < t_min) {
            t_min =  km_ptr->target_pos[i];
        }
        if ( km_ptr -> target_pos[i] > t_max) {
            t_max =  km_ptr->target_pos[i];
        }
        d = (long int) km_ptr->query_pos[i] - (long int) km_ptr->target_pos[i];
        if ( d < d_min ) {
            d_min = d;
        }
        if ( d > d_max ) {
            d_max = d;
        }
    }

    //printf("%lu %ld %ld\n" , km_ptr->count, d_min, d_max);
    if (km_ptr->count == 0)
       d_max = d_min = 0;

    d_count = (seq_coor_t *)calloc( (d_max - d_min)/bin_size + 1, sizeof(seq_coor_t) );
    q_coor =  (seq_coor_t *)calloc( km_ptr->count, sizeof(seq_coor_t) );
    t_coor =  (seq_coor_t *)calloc( km_ptr->count, sizeof(seq_coor_t) );

    for (i = 0; i <  km_ptr->count; i++ ) {
        d = (long int) (km_ptr->query_pos[i]) - (long int) (km_ptr->target_pos[i]);
        d_count[ (d - d_min)/ (long int) bin_size ] += 1;
        q_coor[i] = INT_MAX;
        t_coor[i] = INT_MAX;
    }

    j = 0;
    max_k_mer_count = 0;
    max_k_mer_bin = INT_MAX;
    for (i = 0; i <  km_ptr->count; i++ ) {
        d = (long int) (km_ptr->query_pos[i]) - (long int) (km_ptr->target_pos[i]);
        if ( d_count[ (d - d_min)/ (long int) bin_size ] > max_k_mer_count) {
            max_k_mer_count =  d_count[ (d - d_min)/ (long int) bin_size ];
            max_k_mer_bin = (d - d_min)/ (long int) bin_size;
        }
    }
    //printf("k_mer: %lu %lu\n" , max_k_mer_count, max_k_mer_bin);

    if ( max_k_mer_bin != INT_MAX && max_k_mer_count > count_th ) {
        for (i = 0; i <  km_ptr->count; i++ ) {
            d = (long int) (km_ptr->query_pos[i]) - (long int) (km_ptr->target_pos[i]);
            if ( abs( ( (d - d_min)/ (long int) bin_size ) - max_k_mer_bin ) > 5 ) {
                continue;
            }
            if (d_count[ (d - d_min)/ (long int) bin_size ] > count_th) {
                q_coor[j] = km_ptr->query_pos[i];
                t_coor[j] = km_ptr->target_pos[i];
                //printf("d_count: %lu %lu\n" ,i, d_count[(d - d_min)/ (long int) bin_size]);
                //printf("coor: %lu %lu\n" , q_coor[j], t_coor[j]);
                j ++;
            }
        }
    }

    if (j > 1) {
        arange->s1 = q_coor[0];
        arange->e1 = q_coor[0];
        arange->s2 = t_coor[0];
        arange->e2 = t_coor[0];
        arange->score = 0;

        max_score = 0;
        cur_score = 0;
        cur_start = 0;

        for (i = 1; i < j; i++) {
            cur_score += 32 - (q_coor[i] - q_coor[i-1]);
            //printf("deltaD, %lu %ld\n", q_coor[i] - q_coor[i-1], cur_score);
            if (cur_score < 0) {
                cur_score = 0;
                cur_start = i;
            } else if (cur_score > max_score) {
                arange->s1 = q_coor[cur_start];
                arange->s2 = t_coor[cur_start];
                arange->e1 = q_coor[i];
                arange->e2 = t_coor[i];
                max_score = cur_score;
                arange->score = max_score;
                //printf("%lu %lu %lu %lu\n", arange.s1, arange.e1, arange.s2, arange.e2);
            }
        }

    } else {
        arange->s1 = 0;
        arange->e1 = 0;
        arange->s2 = 0;
        arange->e2 = 0;
        arange->score = 0;
    }

    // printf("free\n");

    free(d_count);
    free(q_coor);
    free(t_coor);
    return arange;
}

aln_range* find_best_aln_range2(kmer_match * km_ptr,
                                seq_coor_t K,
                                seq_coor_t bin_width,
                                seq_coor_t count_th) {

    seq_coor_t * d_coor;
    seq_coor_t * hit_score;
    seq_coor_t * hit_count;
    seq_coor_t * last_hit;
    seq_coor_t max_q, max_t;
    seq_coor_t s, e, max_s, max_e, max_span, d_s, d_e, delta, d_len;
    seq_coor_t px, py, cx, cy;
    seq_coor_t max_hit_idx;
    seq_coor_t max_hit_score, max_hit_count;
    seq_coor_t i, j;
    seq_coor_t candidate_idx, max_d, d;

    aln_range * arange;

    arange = (aln_range *)calloc(1 , sizeof(aln_range));

    d_coor = (seq_coor_t *)calloc( km_ptr->count, sizeof(seq_coor_t) );

    max_q = -1;
    max_t = -1;

    for (i = 0; i <  km_ptr->count; i++ ) {
        d_coor[i] = km_ptr->query_pos[i] - km_ptr->target_pos[i];
        max_q = max_q > km_ptr->query_pos[i] ? max_q : km_ptr->query_pos[i];
        max_t = max_t > km_ptr->target_pos[i] ? max_q : km_ptr->target_pos[i];

    }

    qsort(d_coor, km_ptr->count, sizeof(seq_coor_t), compare_seq_coor);


    s = 0;
    e = 0;
    max_s = -1;
    max_e = -1;
    max_span = -1;
    delta = (long int) ( 0.05 * ( max_q + max_t ) );
    d_len =  km_ptr->count;
    d_s = -1;
    d_e = -1;
    while (1) {
        d_s = d_coor[s];
        d_e = d_coor[e];
        while (d_e < d_s + delta && e < d_len-1) {
            e += 1;
            d_e = d_coor[e];
        }
        if ( max_span == -1 || e - s > max_span ) {
            max_span = e - s;
            max_s = s;
            max_e = e;
        }
        s += 1;
        if (s == d_len || e == d_len) {
            break;
        }
    }

    if (max_s == -1 || max_e == -1 || max_e - max_s < 32) {
        arange->s1 = 0;
        arange->e1 = 0;
        arange->s2 = 0;
        arange->e2 = 0;
        arange->score = 0;
        free(d_coor);
        return arange;
    }

    last_hit = (seq_coor_t *)calloc( km_ptr->count, sizeof(seq_coor_t) );
    hit_score = (seq_coor_t *)calloc( km_ptr->count, sizeof(seq_coor_t) );
    hit_count = (seq_coor_t *)calloc( km_ptr->count, sizeof(seq_coor_t) );

    for (i = 0; i <  km_ptr->count; i++ ) {
        last_hit[i] = -1;
        hit_score[i] = 0;
        hit_count[i] = 0;
    }
    max_hit_idx = -1;
    max_hit_score = 0;
    for (i = 0; i < km_ptr->count; i ++)  {
        cx = km_ptr->query_pos[i];
        cy = km_ptr->target_pos[i];
        d = cx - cy;
        if ( d < d_coor[max_s] || d > d_coor[max_e] ) continue;

        j = i - 1;
        candidate_idx = -1;
        max_d = 65535;
        while (1) {
            if ( j < 0 ) break;
            px = km_ptr->query_pos[j];
            py = km_ptr->target_pos[j];
            d = px - py;
            if ( d < d_coor[max_s] || d > d_coor[max_e] ) {
                j--;
                continue;
            }
            if (cx - px > 320) break; //the number here controling how big alignment gap to be considered
            if (cy > py && cx - px + cy - py < max_d && cy - py <= 320 ) {
                max_d = cx - px + cy - py;
                candidate_idx = j;
            }
            j--;
        }
        if (candidate_idx != -1) {
            last_hit[i] = candidate_idx;
            hit_score[i] = hit_score[candidate_idx] + (64 - max_d);
            hit_count[i] = hit_count[candidate_idx] + 1;
            if (hit_score[i] < 0) {
                hit_score[i] = 0;
                hit_count[i] = 0;
            }
        } else {
            hit_score[i] = 0;
            hit_count[i] = 0;
        }
        if (hit_score[i] > max_hit_score) {
            max_hit_score = hit_score[i];
            max_hit_count = hit_count[i];
            max_hit_idx = i;
        }

    }
    if (max_hit_idx == -1) {
        arange->s1 = 0;
        arange->e1 = 0;
        arange->s2 = 0;
        arange->e2 = 0;
        arange->score = 0;
        free(d_coor);
        free(last_hit);
        free(hit_score);
        free(hit_count);
        return arange;
    }

    arange->score = max_hit_count + 1;
    arange->e1 = km_ptr->query_pos[max_hit_idx];
    arange->e2 = km_ptr->target_pos[max_hit_idx];
    i = max_hit_idx;
    while (last_hit[i] != -1) {
        i = last_hit[i];
    }
    arange->s1 = km_ptr->query_pos[i];
    arange->s2 = km_ptr->target_pos[i];

    free(d_coor);
    free(last_hit);
    free(hit_score);
    free(hit_count);
    return arange;
}

void free_aln_range( aln_range * arange) {
    free(arange);
}
}
