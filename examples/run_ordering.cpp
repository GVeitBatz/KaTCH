/*
 * katch/examples/run_ordering.cpp
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

#include "io/edge_info.h"
#include "io/tpgr_format.h"
#include "io/ddsg_format.h"
#include "prepro/ordering.h"
#include "util/misc.h"

int main(int argc, char** argv)
{
    using Graph = katch::Ordering::Graph;
    using EdgeInfo = katch::EdgeInfo<katch::Ordering::TTF>;

    const char* binary_name = argv[0];

    if ( argc != 4 )
    {
        std::cerr
            << std::endl
            << "USAGE: " << binary_name
            << " <.tpgr input file> <.tch output file> <n threads>" << std::endl
            << std::endl;

        return EXIT_FAILURE;
    }

    const std::string tpgr_input_file_name(argv[1]);
    const std::string tch_output_file_name(argv[2]);
    const int n_threads = std::stoi(std::string(argv[3]));

    std::vector<EdgeInfo> edge_list = katch::tpgr_format::read_edges<EdgeInfo>(tpgr_input_file_name);

    if ( edge_list.empty() )
    {
        KATCH_ERROR("Empty graph.\n");
        return EXIT_FAILURE;
    }

    Graph graph(std::move(edge_list));
    KATCH_STATUS("Graph has " << graph.get_n_nodes() << " nodes, " << graph.get_n_edges() << " edges.\n");

    katch::Ordering ordering(std::move(graph), tch_output_file_name);

    auto t1 = katch::util::time_stamp();
    ordering.order_and_construct(n_threads);
    auto t2 = katch::util::time_stamp();

    KATCH_STATUS("Done.\n");

    return EXIT_SUCCESS;
}
