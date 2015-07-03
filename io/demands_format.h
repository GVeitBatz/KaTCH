/*
 * katch/io/btch_format.h
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

#ifndef KATCH_DEMANDS_FORMAT_H_
#define KATCH_DEMANDS_FORMAT_H_

#include <fstream>
#include <string>
#include <limits>
#include <vector>

#include "datastr/graph/basic.h"
#include "io/binary.h"
#include "util/misc.h"

namespace katch
{

namespace demands_format
{

struct Demand
{
    NodeIterator _start;
    NodeIterator _destination;
    double _t_dep;
    double _t_arr;
    uint32_t _rank;
};

std::vector<Demand> read(const std::string& input_file_name)
{
    std::vector<Demand> demands;

    if (input_file_name == "")
    {
        KATCH_ERROR("Empty input file name given.\n");
        return demands;
    }

    KATCH_STATUS("Reading DEMANDS file '" << input_file_name << "'...");

    std::ifstream input_demands_file(input_file_name);

    if ( ! input_demands_file.is_open() )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Unable to open file.\n");
        return demands;
    }

    char* buffer = new char[10];
    input_demands_file.read(buffer, 9);
    if ( std::string(buffer, 9) != "demands\r\n" )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Not a valid DEMANDS file.\n");
        return demands;
    }

    uint32_t n_demands = read_uint32(input_demands_file);

    for ( size_t i = 0 ; i != n_demands ; ++i )
    {
        uint32_t start = read_uint32(input_demands_file);
        uint32_t destination = read_uint32(input_demands_file);
        double t_dep = read_double(input_demands_file);
        double t_arr = read_double(input_demands_file);
        uint32_t rank = read_uint32(input_demands_file);

        demands.push_back( Demand{ start, destination, t_dep, t_arr, rank } );
    }

    uint32_t terminator = read_uint32(input_demands_file);

    if ( terminator != 0x07162534 )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Not a valid DEMANDS file.\n");
        demands.clear();
        return demands;
    }

    input_demands_file.close();
    KATCH_CONTINUE_STATUS(" OK\n");

    return demands;
}


}
}

#endif /* KATCH_DEMANDS_FORMAT_H_ */
