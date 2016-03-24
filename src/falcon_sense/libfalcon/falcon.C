
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
 *       Filename:  fastcon.c
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

#include "falcon.H"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

namespace FConsensus {

typedef struct {
    seq_coor_t t_pos;
    uint8 delta;
    char q_base;
    seq_coor_t p_t_pos;   // the tag position of the previous base
    uint8 p_delta; // the tag delta of the previous base
    char p_q_base;        // the previous base
    uint32 q_id;
} align_tag_t;

typedef struct {
    seq_coor_t len;
    align_tag_t * align_tags;
} align_tags_t;


typedef struct {
    uint16 size;
    uint16 n_link;
    seq_coor_t * p_t_pos;   // the tag position of the previous base
    uint8 * p_delta; // the tag delta of the previous base
    char * p_q_base;        // the previous base
    uint16 * link_count;
    uint16 count;
    seq_coor_t best_p_t_pos;
    uint8 best_p_delta;
    uint8 best_p_q_base; // encoded base
    double score;
} align_tag_col_t;

typedef struct {
    align_tag_col_t * base;
} msa_base_group_t;

typedef struct {
    uint8 size;
    uint8 max_delta;
    msa_base_group_t * delta;
} msa_delta_group_t;

typedef msa_delta_group_t * msa_pos_t;

align_tags_t * get_align_tags( char * aln_q_seq,
                               char * aln_t_seq,
                               seq_coor_t aln_seq_len,
                               aln_range * range,
                               uint32 q_id,
                               seq_coor_t t_offset) {
    char p_q_base;
    align_tags_t * tags;
    seq_coor_t i, j, jj, k, p_j, p_jj;

    tags = (align_tags_t *)calloc( 1, sizeof(align_tags_t) );
    tags->len = aln_seq_len;
    tags->align_tags = (align_tag_t *)calloc( aln_seq_len + 1, sizeof(align_tag_t) );
    i = range->s1 - 1;
    j = range->s2 - 1;
    jj = 0;
    p_j = -1;
    p_jj = 0;
    p_q_base = '.';

    for (k = 0; k < aln_seq_len; k++) {
        if (aln_q_seq[k] != '-') {
            i ++;
            jj ++;
        }
        if (aln_t_seq[k] != '-') {
            j ++;
            jj = 0;
        }
        //printf("t %d %d %d %c %c\n", q_id, j, jj, aln_t_seq[k], aln_q_seq[k]);


        if ( j + t_offset >= 0 && jj < uint8MAX && p_jj < uint8MAX) {
            (tags->align_tags[k]).t_pos = j + t_offset;
            (tags->align_tags[k]).delta = jj;
            (tags->align_tags[k]).p_t_pos = p_j + t_offset;
            (tags->align_tags[k]).p_delta = p_jj;
            (tags->align_tags[k]).p_q_base = p_q_base;
            (tags->align_tags[k]).q_base = aln_q_seq[k];
            (tags->align_tags[k]).q_id = q_id;

            p_j = j;
            p_jj = jj;
            p_q_base = aln_q_seq[k];
        }
    }
    // sentinal at the end
    //k = aln_seq_len;
    tags->len = k;
    (tags->align_tags[k]).t_pos = uint32MAX;
    (tags->align_tags[k]).delta = uint8MAX;
    (tags->align_tags[k]).q_base = '.';
    (tags->align_tags[k]).q_id = uint32MAX;
    return tags;
}

void free_align_tags( align_tags_t * tags) {
    free( tags->align_tags );
    free( tags );
}


void allocate_aln_col( align_tag_col_t * col) {
    col->p_t_pos = ( seq_coor_t * ) calloc(col->size, sizeof( seq_coor_t ));
    col->p_delta = ( uint8 * ) calloc(col->size, sizeof( uint8 ));
    col->p_q_base = ( char * )calloc(col->size, sizeof( char ));
    col->link_count = ( uint16 * ) calloc(col->size, sizeof( uint16 ));
}

void realloc_aln_col( align_tag_col_t * col ) {
    col->p_t_pos = (seq_coor_t *) realloc( col->p_t_pos, (col->size) * sizeof( seq_coor_t ));
    col->p_delta = ( uint8 *)  realloc( col->p_delta, (col->size) * sizeof( uint8 ));
    col->p_q_base = (char *) realloc( col->p_q_base, (col->size) * sizeof( char ));
    col->link_count = ( uint16 *) realloc( col->link_count, (col->size) * sizeof( uint16 ));
}

void free_aln_col( align_tag_col_t * col) {
    free(col->p_t_pos);
    free(col->p_delta);
    free(col->p_q_base);
    free(col->link_count);
}


void allocate_delta_group( msa_delta_group_t * g) {
    int i,j;
    g->max_delta = 0;
    g->delta = (msa_base_group_t *) calloc( g->size, sizeof(msa_base_group_t));
    for (i = 0; i< g->size; i++) {
        g->delta[i].base = ( align_tag_col_t * ) calloc( 5, sizeof(align_tag_col_t ) );
        for (j = 0; j < 5; j++ ) {
             g->delta[i].base[j].size = 8;
             allocate_aln_col(&(g->delta[i].base[j]));
        }
    }
}

void realloc_delta_group( msa_delta_group_t * g, uint16 new_size ) {
    int i, j, bs, es;
    bs = g->size;
    es = new_size;
    g->delta = (msa_base_group_t *) realloc(g->delta, new_size * sizeof(msa_base_group_t));
    for (i=bs; i < es; i++) {
        g->delta[i].base = ( align_tag_col_t *) calloc( 5, sizeof(align_tag_col_t ) );
        for (j = 0; j < 5; j++ ) {
             g->delta[i].base[j].size = 8;
             allocate_aln_col(&(g->delta[i].base[j]));
        }
    }
    g->size = new_size;
}

void free_delta_group( msa_delta_group_t * g) {
    //manything to do here
    int i, j;
    for (i = 0; i < g->size; i++) {
        for (j = 0; j < 5; j++) {
            free_aln_col( &(g->delta[i].base[j]) );
        }
        free(g->delta[i].base);
    }
    free(g->delta);
}

void update_col( align_tag_col_t * col, seq_coor_t p_t_pos, uint8 p_delta, char p_q_base) {
    int updated = 0;
    int kk;
    col->count += 1;
    for (kk = 0; kk < col->n_link; kk++) {
        if ( p_t_pos == col->p_t_pos[kk] &&
             p_delta == col->p_delta[kk] &&
             p_q_base == col->p_q_base[kk] ) {
            col->link_count[kk] ++;
            updated = 1;
            break;
        }
    }
    if (updated == 0) {
        if (col->n_link + 1 > col->size) {
            if (col->size < (uint16MAX >> 1)-1) {
                col->size *= 2;
            } else {
                col->size += 256;
            }
            assert( col->size < uint16MAX-1 );
            realloc_aln_col(col);
        }
        kk = col->n_link;

        col->p_t_pos[kk] = p_t_pos;
        col->p_delta[kk] = p_delta;
        col->p_q_base[kk] = p_q_base;
        col->link_count[kk] = 1;
        col->n_link++;
    }
}


msa_pos_t * get_msa_working_sapce(uint32 max_t_len) {
    msa_pos_t * msa_array;
    uint32 i;
    msa_array = (msa_pos_t *)calloc(max_t_len, sizeof(msa_pos_t));
    for (i = 0; i < max_t_len; i++) {
        msa_array[i] = (msa_delta_group_t *)calloc(1, sizeof(msa_delta_group_t));
        msa_array[i]->size = 8;
        allocate_delta_group(msa_array[i]);
    }
    return msa_array;
}

void clean_msa_working_space( msa_pos_t * msa_array, uint32 max_t_len) {
    uint32 i,j,k;
    align_tag_col_t * col;
    for (i = 0; i < max_t_len; i++) {
        for (j =0; j < msa_array[i]->max_delta + 1; j++) {
            for (k = 0; k < 5; k++ ) {
                col = msa_array[i]->delta[j].base + k;
                /*
                for (c =0; c < col->size; c++) {
                    col->p_t_pos[c] = 0;
                    col->p_delta[c] = 0;
                    col->p_q_base[c] = 0;
                    col->link_count[c] =0;
                }
                */
                col->n_link = 0;
                col->count = 0;
                col->best_p_t_pos = 0;
                col->best_p_delta = 0;
                col->best_p_q_base = 0;
                col->score = 0;
            }
        }
        msa_array[i]->max_delta = 0;
    }
}

#define STATIC_ALLOCATE
//#undef STATIC_ALLOCATE

consensus_data * get_cns_from_align_tags( align_tags_t ** tag_seqs,
                                          uint32 n_tag_seqs,
                                          uint32 t_len,
                                          uint32 min_cov, uint32 max_len ) {

    seq_coor_t i,j;
    seq_coor_t t_pos = 0;
    seq_coor_t t_count = 0;
    uint32 * coverage;
    uint32 * local_nbase;

    consensus_data * consensus;
    align_tag_t * c_tag;
    static msa_pos_t * msa_array = NULL;

    // figure out true t_len and compact, we might have blank spaces for unaligned sequences
    for (i = 0; i < n_tag_seqs; i++)
       if (tag_seqs[i] != NULL)
          tag_seqs[t_count++] = tag_seqs[i];
    // null out the remainder
    for (i = t_count; i < n_tag_seqs; i++)
       tag_seqs[i] = NULL;
    n_tag_seqs = t_count;

    if (n_tag_seqs == 0) {
        // allocate an empty consensus sequence
        consensus = (consensus_data *)calloc( 1, sizeof(consensus_data) );
        consensus->sequence = (char *)calloc( 1, sizeof(char) );
        consensus->eqv = (int32 *)calloc( 1, sizeof(int32) );
        return consensus;
    }

    coverage = (uint32 *)calloc( t_len, sizeof(uint32) );
    local_nbase = (uint32 *)calloc( t_len, sizeof(uint32) );

#ifndef STATIC_ALLOCATE

    msa_array = (msa_pos_t *)calloc(t_len, sizeof(msa_pos_t));

    for (i = 0; i < t_len; i++) {
        msa_array[i] = (msa_delta_group_t *)calloc(1, sizeof(msa_delta_group_t));
        msa_array[i]->size = 8;
        allocate_delta_group(msa_array[i]);
    }

#endif

#ifdef STATIC_ALLOCATE

    if ( msa_array == NULL) {
        msa_array = get_msa_working_sapce( max_len );
    }

    assert(t_len < max_len);

#endif


    // loop through every alignment
    //printf("XX %d\n", n_tag_seqs);
    for (i = 0; i < n_tag_seqs; i++) {

        // for each alignment position, insert the alignment tag to msa_array
        for (j = 0; j < tag_seqs[i]->len; j++) {
            c_tag = tag_seqs[i]->align_tags + j;
            uint32 delta;
            delta = c_tag->delta;
            if (delta == 0) {
                t_pos = c_tag->t_pos;
                coverage[ t_pos ] ++;
            }
            // Assume t_pos was set on earlier iteration.
            // (Otherwise, use its initial value, which might be an error. ~cd)
            if (delta > msa_array[t_pos]->max_delta) {
                msa_array[t_pos]->max_delta = delta;
                if (msa_array[t_pos]->max_delta + 4 > msa_array[t_pos]->size ) {
                    realloc_delta_group(msa_array[t_pos], msa_array[t_pos]->max_delta + 8);
                }
            }

            uint32 base = -1;
            switch (c_tag->q_base) {
                case 'A': base = 0; break;
                case 'C': base = 1; break;
                case 'G': base = 2; break;
                case 'T': base = 3; break;
                case '-': base = 4; break;
                default : base = 4; break;
            }
            // Note: On bad input, base may be -1.
            update_col( &(msa_array[t_pos]->delta[delta].base[base]), c_tag->p_t_pos, c_tag->p_delta, c_tag->p_q_base);
            local_nbase[ t_pos ] ++;
        }
    }

    // propogate score throught the alignment links, setup backtracking information
    align_tag_col_t * g_best_aln_col = 0;
    uint32 g_best_ck = 0;
    seq_coor_t g_best_t_pos = 0;
    {
        int kk;
        int ck;
        // char base;
        int best_i;
        int best_j;
        int best_b;
        int best_ck = -1;
        double score;
        double best_score;
        double g_best_score;
        // char best_mark;

        align_tag_col_t * aln_col;

        g_best_score = -1;

        for (i = 0; i < t_len; i++) {  //loop through every template base
            //printf("max delta: %d %d\n", i, msa_array[i]->max_delta);
            for (j = 0; j <= msa_array[i]->max_delta; j++) { // loop through every delta position
                for (kk = 0; kk < 5; kk++) {  // loop through diff bases of the same delta posiiton
                    /*
                    switch (kk) {
                        case 0: base = 'A'; break;
                        case 1: base = 'C'; break;
                        case 2: base = 'G'; break;
                        case 3: base = 'T'; break;
                        case 4: base = '-'; break;
                    }
                    */
                    aln_col = msa_array[i]->delta[j].base + kk;
                    best_score = -1;
                    best_i = -1;
                    best_j = -1;
                    best_b = -1;

                    for (ck = 0; ck < aln_col->n_link; ck++) { // loop through differnt link to previous column
                        int pi;
                        int pj;
                        int pkk;
                        pi = aln_col->p_t_pos[ck];
                        pj = aln_col->p_delta[ck];
                        switch (aln_col->p_q_base[ck]) {
                            case 'A': pkk = 0; break;
                            case 'C': pkk = 1; break;
                            case 'G': pkk = 2; break;
                            case 'T': pkk = 3; break;
                            case '-': pkk = 4; break;
                            default : pkk = 4; break;
                        }

                        if (aln_col->p_t_pos[ck] == -1) {
                            score =  (double) aln_col->link_count[ck] - (double) coverage[i] * 0.5;
                        } else if (pj > msa_array[pi]->max_delta) {
                            score =  (double) aln_col->link_count[ck] - (double) coverage[i] * 0.5;
                        } else {
                            score = msa_array[pi]->delta[pj].base[pkk].score +
                                    (double) aln_col->link_count[ck] - (double) coverage[i] * 0.5;
                        }
                        // best_mark = ' ';
                        if (score > best_score) {
                            best_score = score;
                            aln_col->best_p_t_pos = best_i = pi;
                            aln_col->best_p_delta = best_j = pj;
                            aln_col->best_p_q_base = best_b = pkk;
                            best_ck = ck;
                            // best_mark = '*';
                        }
                        /*
                        printf("X %d %d %d %c %d %d %d %c %d %lf %c\n", coverage[i], i, j, base, aln_col->count,
                                                              aln_col->p_t_pos[ck],
                                                              aln_col->p_delta[ck],
                                                              aln_col->p_q_base[ck],
                                                              aln_col->link_count[ck],
                                                              score, best_mark);
                        */
                    }
                    aln_col->score = best_score;
                    if (best_score > g_best_score) {
                        g_best_score = best_score;
                        g_best_aln_col = aln_col;
                        g_best_ck = best_ck;
                        g_best_t_pos = i;
                        //printf("GB %d %d %d %d\n", i, j, ck, g_best_aln_col);
                    }
                }
            }
        }
        assert(g_best_score != -1);
    }

    // reconstruct the sequences
    uint32 index;
    char bb = '$';
    int ck;
    char * cns_str;
    int * eqv;
    double score0;

    consensus = (consensus_data *)calloc( 1, sizeof(consensus_data) );
    consensus->sequence = (char *)calloc( t_len * 2 + 1, sizeof(char) );
    consensus->eqv = (int32 *)calloc( t_len * 2 + 1, sizeof(int32) );
    cns_str = consensus->sequence;
    eqv =  consensus->eqv;

    index = 0;
    ck = g_best_ck;
    i = g_best_t_pos;

    while (1) {
        if (coverage[i] > min_cov) {
            switch (ck) {
                case 0: bb = 'A'; break;
                case 1: bb = 'C'; break;
                case 2: bb = 'G'; break;
                case 3: bb = 'T'; break;
                case 4: bb = '-'; break;
            }
        } else {
            switch (ck) {
                case 0: bb = 'a'; break;
                case 1: bb = 'c'; break;
                case 2: bb = 'g'; break;
                case 3: bb = 't'; break;
                case 4: bb = '-'; break;
            }
        }
        // Note: On bad input, bb will keep previous value, possibly '$'.

        score0 = g_best_aln_col->score;
        i = g_best_aln_col->best_p_t_pos;
        if (i == -1 || index >= t_len * 2) break;
        j = g_best_aln_col->best_p_delta;
        ck = g_best_aln_col->best_p_q_base;
        g_best_aln_col = msa_array[i]->delta[j].base + ck;

        if (bb != '-') {
            cns_str[index] = bb;
            eqv[index] = (int) score0 - (int) g_best_aln_col->score;
            //printf("C %d %d %c %lf %d %d\n", i, index, bb, g_best_aln_col->score, coverage[i], eqv[index] );
            index ++;
        }
    }

    // reverse the sequence
    for (i = 0; i < index/2; i++) {
        cns_str[i] = cns_str[i] ^ cns_str[index-i-1];
        cns_str[index-i-1] = cns_str[i] ^ cns_str[index-i-1];
        cns_str[i] = cns_str[i] ^ cns_str[index-i-1];
        eqv[i] = eqv[i] ^ eqv[index-i-1];
        eqv[index-i-1] = eqv[i] ^ eqv[index-i-1];
        eqv[i] = eqv[i] ^ eqv[index-i-1];
    }

    cns_str[index] = 0;
    //printf("%s\n", cns_str);
#ifndef STATIC_ALLOCATE
    for (i = 0; i < t_len; i++) {
        free_delta_group(msa_array[i]);
        free(msa_array[i]);
    }

    free(msa_array);
#endif

#ifdef STATIC_ALLOCATE
    clean_msa_working_space(msa_array, t_len+1);
#endif

    free(coverage);
    free(local_nbase);
    return consensus;
}

consensus_data * generate_consensus( vector<string> input_seq,
                           uint32 min_cov,
                           uint32 K,
                           double min_idt, uint32 min_len, uint32 max_len) {
    uint32 seq_count;
    kmer_lookup * lk_ptr;
    seq_array sa_ptr;
    seq_addr_array sda_ptr;
    align_tags_t ** tags_list;
    consensus_data * consensus;
    double max_diff;
    max_diff = 1.0 - min_idt;

    seq_count = input_seq.size();
    fflush(stdout);

    tags_list = (align_tags_t **)calloc( seq_count, sizeof(align_tags_t*) );
    lk_ptr = allocate_kmer_lookup( 1 << (K * 2) );
    sa_ptr = allocate_seq( (seq_coor_t) input_seq[0].length() );
    sda_ptr = allocate_seq_addr( (seq_coor_t) input_seq[0].length() );
    add_sequence( 0, K, input_seq[0].c_str(), input_seq[0].length(), sda_ptr, sa_ptr, lk_ptr);

#pragma omp parallel for schedule(dynamic)
    for (uint32 j=0; j < seq_count; j++) {
#define MAX_UNMASKED_LENGTH 500000
#define MAX_KMER_REPEAT     1000
        if (input_seq[j].length() > MAX_UNMASKED_LENGTH) {
            mask_k_mer(1 << (K*2), lk_ptr, MAX_KMER_REPEAT);
        }
        kmer_match *kmer_match_ptr = find_kmer_pos_for_seq(input_seq[j].c_str(), input_seq[j].length(), K, sda_ptr, lk_ptr);
#define INDEL_ALLOWENCE_0 6

        aln_range *arange = find_best_aln_range(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels

        //fprintf(stderr, "1:read %d %ld %ld %ld %ld\n", j, arange->s1, arange->e1, arange->s2, arange->e2);

        //arange = find_best_aln_range2(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels

        //printf("2:%ld %ld %ld %ld\n\n", arange->s1, arange->e1, arange->s2, arange->e2);

#define INDEL_ALLOWENCE_1 0.10
        if (arange->e1 - arange->s1 < 100 || arange->e2 - arange->s2 < 100 ||
            abs( (arange->e1 - arange->s1 ) - (arange->e2 - arange->s2) ) >
                   (int) (0.5 * INDEL_ALLOWENCE_1 * (arange->e1 - arange->s1 + arange->e2 - arange->s2))) {
            free_kmer_match( kmer_match_ptr);
            free_aln_range(arange);
            continue;
        }


#define INDEL_ALLOWENCE_2 150
        NDalignment::NDalignResult aln;
        align(input_seq[j].c_str()+arange->s1, arange->e1 - arange->s1 ,
                    input_seq[0].c_str()+arange->s2, arange->e2 - arange->s2 ,
                    INDEL_ALLOWENCE_2, 1, aln);
        if (aln._size > min_len && ((double) aln._dist / (double) aln._size) < max_diff) {
            tags_list[j] = get_align_tags( aln._qry_aln_str,
                                                           aln._tgt_aln_str,
                                                           aln._size,
                                                           arange, j,
                                                           0);
           //fprintf(stderr, "Aligned seq %d  to positions %d - %d and %d - %d with %d diffs size %d\n", j, aln._qry_bgn, aln._qry_end, aln._tgt_bgn, aln._tgt_end, aln._dist, aln._size);
        }
        free_aln_range(arange);
        free_kmer_match( kmer_match_ptr);
    }

    consensus = get_cns_from_align_tags( tags_list, seq_count, input_seq[0].length(), min_cov, max_len);
    free_seq_addr_array(sda_ptr);
    free_seq_array(sa_ptr);
    free_kmer_lookup(lk_ptr);
    for (int j=0; j < seq_count; j++)
        if (tags_list[j] != NULL)
           free_align_tags(tags_list[j]);
    free(tags_list);
    return consensus;
}

void free_consensus_data( consensus_data * consensus ){
    free(consensus->sequence);
    free(consensus->eqv);
    free(consensus);
}
}
