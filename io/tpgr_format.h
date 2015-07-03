/*
 * katch/io/tpgr_format.h
 *
 *
 * This file is part of
 *
 * KaTCH -- Karlsruhe Time-Dependent Contraction Hierarchies
 *
 * Copyright (C) 2015
 *
 * Institut fuer Theroretische Informatik,
 * Karlsruher Institut fuer Technology (KIT),
 * 76131 Karlsruhe, Germany
 *
 * Author: Gernot Veit Batz
 *
 *
 *
 * KaTCH is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaTCH is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with Contraction Hierarchies; see the file COPYING;
 * if not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef KATCH_TPGR_FORMAT_H_
#define KATCH_TPGR_FORMAT_H_

#include <fstream>
#include <vector>
#include <string>
#include <utility>

#include "io/edge_info.h"
#include "datastr/base/point.h"
#include "datastr/base/ttf_wrapper.h"
#include "util/misc.h"

namespace katch
{

namespace tpgr_format
{

template <typename EdgeInfo>
std::vector<EdgeInfo> read_edges(const std::string& input_file_name)
{
    using TTF = typename EdgeInfo::TTF;

    std::vector<EdgeInfo> result;

    size_t check_n_edges = 0;
    size_t check_n_points = 0;

    if (input_file_name == "")
    {
        KATCH_ERROR("Empty input file name given.\n");
        return result;
    }

    KATCH_STATUS("Reading TPGR file '" << input_file_name << "'...");

    std::ifstream input_tpgr_file(input_file_name);

    if ( ! input_tpgr_file.is_open() )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Unable to open file '" << input_file_name << "'\n");
        return result;
    }

    unsigned int n_nodes;
    input_tpgr_file >> n_nodes;

    unsigned int n_edges;
    input_tpgr_file >> n_edges;

    unsigned int n_points;
    input_tpgr_file >> n_points;

    unsigned int period;
    input_tpgr_file >> period;

    unsigned int edge_counter = 0;
    while ( ! input_tpgr_file.eof() && edge_counter < n_edges )
    {
        ++check_n_edges;
        ++edge_counter;

        unsigned int src;
        input_tpgr_file >> src;

        unsigned int tgt;
        input_tpgr_file >> tgt;

        unsigned int sample_size;
        input_tpgr_file >> sample_size;

        if ( sample_size == 0 )
        {
            KATCH_CONTINUE_STATUS(" ABORT\n");
            KATCH_ERROR("TPGR file corrupted: Time function has zero points.\n");
            result.clear();
            return result;
        }

        std::vector<Point> sample_points;
        for ( size_t i = 0; i < sample_size ; i++ )
        {
            ++check_n_points;

            double x;
            input_tpgr_file >> x;

            double y;
            input_tpgr_file >> y;

            if ( x < 0 || x >= period )
            {
                KATCH_CONTINUE_STATUS(" ABORT\n");
                KATCH_ERROR("TPGR file corrupted: x-value not in [0,period).\n");
                result.clear();
                return result;
            }

            if ( y < 0 )
            {
                KATCH_CONTINUE_STATUS(" ABORT\n");
                KATCH_ERROR("TPGR file corrupted: y-value smaller than 0.\n");
                result.clear();
                return result;
            }

            if ( ! sample_points.empty() && x < sample_points.back().x )
            {
                KATCH_CONTINUE_STATUS(" ABORT\n");
                KATCH_ERROR("TPGR file corrupted: x-value smaller than the one before.\n");
                result.clear();
                return result;
            }

            sample_points.push_back(Point(x, y));
        }

        assert( sample_points.size() > 0 );

        EdgeInfo edge_info;
        edge_info.set_source(src);
        edge_info.set_target(tgt);
        edge_info.set_forward(true);
        edge_info.set_backward(false);
        edge_info.set_ttf( TTF(sample_points.begin(), sample_points.end()) );

        result.push_back(std::move(edge_info));
    }

    if ( check_n_edges != n_edges )
    {

        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("TPGR file corrupted: wrong number of edges.\n");
        result.clear();
        return result;
   }

    if ( check_n_points != n_points )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("TPGR file corrupted: wrong number of bend points\n");
        result.clear();
        return result;
    }

    input_tpgr_file.close();
    KATCH_CONTINUE_STATUS(" OK\n");

    return result;
}

}

}

#endif /* KATCH_TPGR_FORMAT_H_ */
