/*
 * katch/examples/run_tch_ea_query.cpp
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

#include <iostream>
#include <stdlib.h>
#include <string>

#include "util/misc.h"
#include "io/btch_format.h"
#include "io/demands_format.h"
#include "io/edge_info.h"
#include "query/tch_ea_query.h"

int main(int argc, char** argv)
{
    using Graph = katch::TchEaQuery::Graph;
    using EdgeInfo = katch::EdgeInfo<Graph::TTF>;

    const char* binary_name = argv[0];

    if ( argc != 3 )
    {
        std::cerr
            << std::endl
            << "USAGE: " << binary_name
            << " <.btch input file> <.demands input file>" << std::endl
            << std::endl;

        return EXIT_FAILURE;
    }

    const std::string btch_input_file_name(argv[1]);
    const std::string demands_input_file_name(argv[2]);

    std::vector<EdgeInfo> edge_list = katch::btch_format::read_edges<EdgeInfo>(btch_input_file_name);

    Graph graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");

    std::vector<katch::demands_format::Demand> demands = katch::demands_format::read(demands_input_file_name);
    katch::TchEaQuery query(std::move(graph));

    KATCH_STATUS("Running " << demands.size() << " queries...");

    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    double total_abs_error = 0.0;
    double total_rel_error = 0.0;

    double max_abs_error_expand = 0.0;
    double max_rel_error_expand = 0.0;
    double total_abs_error_expand = 0.0;
    double total_rel_error_expand = 0.0;

    double total_route_extraction_and_expansion_time = 0.0;

    auto t1 = katch::util::time_stamp();

    for ( const auto& demand : demands )
    {
        // obtain travel time from start to destination
        const double t_arr = query.one_to_one(demand._start, demand._destination, demand._t_dep);

        // compute absolute and relative errors of obtained travel time
        const double abs_error = fabs(t_arr - demand._t_arr);
        max_abs_error = std::max(max_abs_error, abs_error);
        total_abs_error += abs_error;

        const double rel_error = abs_error / demand._t_arr;
        max_rel_error = std::max(max_rel_error, rel_error);
        total_rel_error += rel_error;

        auto t_expand_begin = katch::util::time_stamp();

        // obtain up-down-path from start to destination (requires katch::TchEaQuery::one_to_one called first)
        const katch::TchEaQuery::Path up_down_path = query.get_up_down_path(demand._destination);

        // expand up-down-path to obtain path in the original road network
        const katch::TchEaQuery::Path original_path = query.expand_path(up_down_path);

        auto t_expand_end = katch::util::time_stamp();
        total_route_extraction_and_expansion_time += katch::util::get_duration_in_seconds(t_expand_begin, t_expand_end);

        // compute absolute and relative errors of the travel time of the expanded path
        const double t_arr_expand = original_path.get_t_arr(original_path.get_n_edges());

        const double abs_error_expand = fabs(t_arr_expand - demand._t_arr);
        max_abs_error_expand = std::max(max_abs_error_expand, abs_error_expand);
        total_abs_error_expand += abs_error_expand;

        const double rel_error_expand = abs_error_expand / demand._t_arr;
        max_rel_error_expand = std::max(max_rel_error_expand, rel_error_expand);
        total_rel_error_expand += rel_error_expand;
    }

    auto t2 = katch::util::time_stamp();

    KATCH_CONTINUE_STATUS(" OK\n");

    KATCH_STATUS("avg. running time (including route extraction) = " << ((katch::util::get_duration_in_seconds(t1, t2) / double(demands.size())) *1000.0) << " msec\n");
    KATCH_STATUS("avg. running time (without route extraction)   = " << (((katch::util::get_duration_in_seconds(t1, t2) - total_route_extraction_and_expansion_time) / double(demands.size())) *1000.0) << " msec\n");
    KATCH_STATUS("\n");
    KATCH_STATUS("Error of computed arrival time\n");
    KATCH_STATUS("------------------------------\n");
    KATCH_STATUS("max. rel. error = " << (max_rel_error * 100.0) << " %\n");
    KATCH_STATUS("avg. rel. error = " << ((total_rel_error / double(demands.size())) * 100.0) << " %\n");

    KATCH_STATUS("\n");
    KATCH_STATUS("Error of arrival time of expanded route\n");
    KATCH_STATUS("---------------------------------------\n");
    KATCH_STATUS("max. rel. error = " << (max_rel_error_expand  * 100.0) << " %\n");
    KATCH_STATUS("avg. rel. error = " << ((total_rel_error_expand / double(demands.size())) * 100.0) << " %\n");
//    KATCH_STATUS("avg. running time = " << ((katch::util::get_duration_in_seconds(t1, t2)  / double(demands.size())) *1000.0) << " msec\n");

    return EXIT_SUCCESS;
}
