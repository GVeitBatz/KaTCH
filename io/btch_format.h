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

#ifndef KATCH_BTCH_FORMAT_H_
#define KATCH_BTCH_FORMAT_H_

#include <algorithm>
#include <fstream>
#include <string>
#include <limits>
#include <vector>

#include "datastr/base/point.h"
#include "datastr/graph/basic.h"
#include "io/binary.h"
#include "util/misc.h"

namespace katch
{

namespace btch_format
{

constexpr uint32_t FILE_SEPARATOR = 0x00772255;
constexpr uint32_t FILE_TERMINATOR = 0xaabbccdd;
constexpr uint32_t UINT32_MAX_VALUE = std::numeric_limits<uint32_t>::max();

class OutputFile
{

private:

    std::ofstream _file_output_stream;
    std::string _file_name;
    uint32_t _period;
    uint32_t _n_nodes;
    uint32_t _current_n_edges;

    void write_uint32(const uint32_t value)
    {
        katch::write_uint32(value, _file_output_stream);
    }

    void write_double(const double value)
    {
        katch::write_double(value, _file_output_stream);
    }

public:

    OutputFile(const std::string& file_name, const uint32_t& n_nodes, const uint32_t& period)
    : _file_output_stream( file_name, std::fstream::binary ),
      _file_name(file_name),
      _period( period ),
      _n_nodes( n_nodes ),
      _current_n_edges( 0 )
    {
        if ( _file_name.empty() )
            KATCH_ERROR("Empty output file name given.\n");

        if ( _file_output_stream.good() )
        {
            // write format tag
            _file_output_stream.write("BTCH\r\n", 6);

            // write version number
            write_uint32(1);

            // write number of nodes
            write_uint32(n_nodes);

            // write a place holder for the number of edges
            write_uint32(UINT32_MAX_VALUE);

            // write period of TTFs
            write_uint32(period);

            // write a place holder for the level info
            for ( uint32_t node_id = 0; node_id < n_nodes; ++node_id )
                write_uint32(UINT32_MAX_VALUE);

            // write separator between level info and edge info
            write_uint32(FILE_SEPARATOR);
        }
        else
        {
            KATCH_ERROR("BTCH output file '" << _file_name << "' not good.\n");
        }
    }

    ~OutputFile()
    {
        if ( _file_output_stream.is_open() )
        {
            KATCH_WARNING("Automatically closing non-closed BTCH output file on destruction.\n");
            close();
        }
    }

    void write_level_info(const std::vector<uint32_t> node_id_to_level)
    {
        if ( ! _file_output_stream.good() )
        {
            KATCH_ERROR("BTCH output file '" << _file_name << "' not good.\n");
            return;
        }

        if ( _n_nodes != node_id_to_level.size() )
        {
            KATCH_ERROR("Node level information has wrong size.\n");
            return;
        }

        if ( util::has_duplicates(node_id_to_level) )
        {
            KATCH_ERROR("Non-unique node levels encountered.\n");
            return;
        }

        auto old_pos = _file_output_stream.tellp();
        _file_output_stream.seekp( 6 + 4 * sizeof(uint32_t) );

        for ( uint32_t u = 0 ; u < _n_nodes ; ++u )
        {
            assert( u < node_id_to_level.size() );
            const uint32_t level = node_id_to_level[u];

            if ( level >= _n_nodes )
            {
                KATCH_ERROR("Level exceeds number of nodes (" << level << " >= " << _n_nodes << ").\n");
                return;
            }

            write_uint32(level);
        }

        _file_output_stream.seekp( old_pos );
    }

    template <typename OrderEdgeDataT>
    void write_edge(const uint32_t source_node_id, const uint32_t target_node_id, const OrderEdgeDataT& edge_data)
    {
        if ( ! _file_output_stream.good() )
        {
            KATCH_ERROR("BTCH output file '" << _file_name << "' not good.\n");
            return;
        }

        write_uint32(source_node_id);
        write_uint32(target_node_id);

        // write shortcut descriptor
        if ( edge_data.never_represents_a_shortcut() )
            write_uint32(0);
        else if ( edge_data.get_n_middle_nodes() == 1 )
        {
            write_uint32(1);
            write_double(0.0);
            write_uint32( edge_data.get_middle_node_descr(0)._middle_node );
        }
        else
        {
            write_uint32( edge_data.get_n_middle_nodes() );

            for ( size_t i = 0; i < edge_data.get_n_middle_nodes() ; ++i )
            {
                write_double( edge_data.get_middle_node_descr(i)._time  );

                if ( edge_data.get_middle_node_descr(i)._middle_node == INVALID_NODE_ITERATOR )
                    write_uint32(UINT32_MAX_VALUE);
                else
                    write_uint32( edge_data.get_middle_node_descr(i)._middle_node );
            }
        }

        // write travel time function
        if ( edge_data.get_ttf().is_constant() )
        {
            write_uint32(1);
            write_double(0.0);
            write_double(edge_data.get_ttf().get_constant_value());
        }
        else
        {
            assert( edge_data.get_ttf().size() > 0 );

            write_uint32( edge_data.get_ttf().size() );
            for ( size_t i = 0 ; i < edge_data.get_ttf().size() ; i++ )
            {
                write_double( edge_data.get_ttf()[i].x );
                write_double( edge_data.get_ttf()[i].y );
            }
        }

        ++_current_n_edges;
    }

    void close()
    {
        if ( ! _file_output_stream.is_open() )
        {
            KATCH_ERROR("Trying to close a non-open BTCH file '" << _file_name << "'.");
            return;
        }

        if ( ! _file_output_stream.good() )
        {
            KATCH_ERROR("BTCH output file '" << _file_name << "' not good.\n");
            return;
        }

        // append file terminator
        write_uint32(FILE_TERMINATOR);

        // store number of edges
        _file_output_stream.seekp( 6 + 2 * sizeof(uint32_t) );
        write_uint32(_current_n_edges);

        _file_output_stream.close();
    }


};


template <typename EdgeInfo>
std::vector<EdgeInfo> read_edges(const std::string& input_file_name, double epsilon = 0.0)
{
    using TTF = typename EdgeInfo::TTF;

    std::vector<EdgeInfo> result;

    if (input_file_name == "")
    {
        KATCH_ERROR("Empty input file name given.\n");
        return result;
    }

    KATCH_STATUS("Reading BTCH file '" << input_file_name << "'...");

    std::ifstream input_btch_file(input_file_name);

    if ( ! input_btch_file.is_open() )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Unable to open file.\n");
        return result;
    }

    char* buffer = new char[8];
    input_btch_file.read(buffer, 6);
    if ( std::string(buffer, 6) != "BTCH\r\n" )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Not a valid BTCH file.\n");
        return result;
    }

    uint32_t version = read_uint32(input_btch_file);

    if ( version != 1 )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("file has wrong version (expected 1).\n");
        return result;
    }

    uint32_t n_nodes = read_uint32(input_btch_file);
    uint32_t n_edges = read_uint32(input_btch_file);
    double period = read_uint32(input_btch_file);

    if ( period != 864000 )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("file has wrong period of TTFs (expected 864000).\n");
        return result;
    }
    
    std::vector<uint32_t> level;
    level.reserve(n_nodes);

    for ( size_t i = 0 ; i != n_nodes ; ++i )
        level.push_back(read_uint32(input_btch_file));

    if ( util::has_duplicates(level.begin(), level.end()) )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Not a valid BTCH file.\n");
        return result;
    }

    uint32_t separator = read_uint32(input_btch_file);

    if ( separator != FILE_SEPARATOR )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Not a valid BTCH file.\n");
        return result;
    }

    for ( size_t i = 0 ; i != n_edges ; ++i )
    {
        uint32_t source = read_uint32(input_btch_file);
        uint32_t target = read_uint32(input_btch_file);
        uint32_t n_middle_nodes = read_uint32(input_btch_file);

        std::vector<MiddleNodeDescr> middle_node_data;
        middle_node_data.reserve(n_middle_nodes);

        for ( size_t j = 0 ; j != n_middle_nodes ; ++j )
        {
            double time = read_double(input_btch_file);
            uint32_t middle_node = read_uint32(input_btch_file);

            NodeIterator middle_node_it =
                    middle_node == UINT32_MAX_VALUE ?
                            INVALID_NODE_ITERATOR :
                            NodeIterator(middle_node);

            middle_node_data.emplace_back(time, middle_node_it);
        }

        uint32_t n_points = read_uint32(input_btch_file);

        std::vector<Point> sample_points;
        for ( size_t j = 0 ; j != n_points ; ++j )
        {
            double x = read_double(input_btch_file);
            double y = read_double(input_btch_file);

            if ( x < 0 || x > period )
            {
                KATCH_CONTINUE_STATUS(" ABORT\n");
                KATCH_ERROR("BTCH file corrupted: x-value not in [0,period).\n");
                result.clear();
                return result;
            }

            if ( y < 0 )
            {
                KATCH_CONTINUE_STATUS(" ABORT\n");
                KATCH_ERROR("BTCH file corrupted: y-value smaller than 0.\n");
                result.clear();
                return result;
            }

            if ( ! sample_points.empty() > 0 && x < sample_points.back().x )
            {
                KATCH_CONTINUE_STATUS(" ABORT\n");
                KATCH_ERROR("BTCH file corrupted: x-value smaller than the one before.\n");
                result.clear();
                return result;
            }

            sample_points.push_back(Point(x, y));
        }

        EdgeInfo edge_info;

        assert( source < level.size() );
        assert( target < level.size() );
        if ( level[source] < level[target] )
        {
            edge_info.set_source(source);
            edge_info.set_target(target);
            edge_info.set_forward(true);
            edge_info.set_backward(false);
        }
        else
        {
            assert( level[source] > level[target] );

            edge_info.set_source(target);
            edge_info.set_target(source);
            edge_info.set_forward(false);
            edge_info.set_backward(true);
        }

        assert( sample_points.size() > 0 );
        edge_info.set_ttf( TTF(sample_points.begin(), sample_points.end()) );
        edge_info.set_middle_node_data( std::move(middle_node_data) );

        result.push_back(std::move(edge_info));
    }

    uint32_t terminator = read_uint32(input_btch_file);

    if ( terminator != FILE_TERMINATOR )
    {
        KATCH_CONTINUE_STATUS(" ABORT\n");
        KATCH_ERROR("Not a valid BTCH file.\n");
        return result;
    }

    input_btch_file.close();
    KATCH_CONTINUE_STATUS(" OK\n");

    return result;
}


}
}

#endif /* KATCH_BTCH_FORMAT_H_ */
