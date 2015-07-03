/*
 * katch/processing/preprocessing.h
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

#ifndef KATCH_PREPROCESSING_H_
#define KATCH_PREPROCESSING_H_

#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <queue>

#include "datastr/base/ttf_wrapper.h"
#include "datastr/graph/dynamic_graph.h"
#include "datastr/graph/order_edge_data.h"
#include "datastr/graph/basic.h"
#include "datastr/base/double.h"
#include "datastr/base/interval.h"
#include "io/btch_format.h"
#include "util/misc.h"

#include "prepro/edge_description.h"
#include "prepro/witness_search.h"

namespace katch
{

class Ordering
{
public:

    using Graph = DynamicGraph<OrderEdgeData>;
    using EdgeData = Graph::EdgeData;
    using TTF = OrderEdgeData::TTF;
    using UndercutDescriptor = TTF::UndercutDescriptor;

    static constexpr int MAXIMUM_THREADS = -1;

private:

    static bool constexpr SIMULATE = true;
    static bool constexpr DO_NOT_SIMULATE = false;

    using ShortcutStatus = WitnessSearch::ShortcutStatus;

    class WitnessCacheEntry
    {
    public:

        ShortcutStatus _shortcut_status : 2;
        std::uint32_t _shortcut_complexity : 30;
        NodeIterator _u;
        NodeIterator _x;
        NodeIterator _v;

        WitnessCacheEntry()
        : _shortcut_status(WitnessSearch::UNDECIDED), _shortcut_complexity(0),
          _u(INVALID_NODE_ITERATOR), _x(INVALID_NODE_ITERATOR), _v(INVALID_NODE_ITERATOR)
        {}

        WitnessCacheEntry
        (
                const ShortcutStatus& shortcut_status,
                const std::uint32_t& shortcut_complexity,
                const NodeIterator& u, const NodeIterator& x, const NodeIterator& v
        )
        : _shortcut_status(shortcut_status), _shortcut_complexity(shortcut_complexity),
          _u(u), _x(x), _v(v)
        {}
    };

    class WitnessCache
    {
    private:
        std::vector<std::vector<WitnessCacheEntry>> _data;

    public:

        WitnessCache(size_t n_nodes)
        : _data(n_nodes)
        {}

        static const WitnessCacheEntry INVALID_ENTRY;

        const WitnessCacheEntry& lookup(const NodeIterator& u, const NodeIterator& x, const NodeIterator& v) const
        {
            for ( const auto& cache_entry : _data[x] )
            {
                assert(x == cache_entry._x);
                if ( cache_entry._u == u && cache_entry._v == v ) return cache_entry;
            }

            return INVALID_ENTRY;
        }

        bool contains(const NodeIterator& u, const NodeIterator& x, const NodeIterator& v) const
        {
            assert( KATCH_EQUIV( lookup(u, x, v)._x == INVALID_NODE_ITERATOR, &lookup(u, x, v) == &INVALID_ENTRY) );
            return lookup(u, x, v)._x != INVALID_NODE_ITERATOR;
        }

        void insert(const WitnessCacheEntry& cache_entry)
        {
            assert( ! contains(cache_entry._u, cache_entry._x, cache_entry._v) );
            _data[cache_entry._x].push_back(cache_entry);
        }

        void remove(NodeIterator x)
        {
            _data[x].clear();
            _data[x].shrink_to_fit();
        }

        void remove(const NodeIterator& u, const NodeIterator& x, const NodeIterator& v)
        {
            for ( int i = 0 ; i < (int) _data[x].size() ; ++i )
            {
                if
                (
                        ( u == INVALID_NODE_ITERATOR || u == _data[x][i]._u ) &&
                        ( v == INVALID_NODE_ITERATOR || v == _data[x][i]._v)
                )
                {
                    _data[x][i] = _data[x].back();
                    _data[x].pop_back();
                    --i;
                }
            }
        }
    };

    class LocalThreadData
    {

    public:

        std::vector<EdgeDescription> _edges_to_insert;
        std::vector<WitnessCacheEntry> _witnesses_to_insert_to_cache;
        WitnessSearch _witness_search;

        LocalThreadData(const LocalThreadData&) = delete;
        LocalThreadData& operator=(const LocalThreadData&) = delete;

        LocalThreadData(LocalThreadData&&) = default;
        LocalThreadData& operator=(LocalThreadData&&) = default;

        LocalThreadData(const Graph* graph)
        : _edges_to_insert(),
          _witnesses_to_insert_to_cache(),
          _witness_search(graph)
        {}
    };

    struct ContractJob
    {
        BwEdgeIterator _e_in;
        NodeIterator _x;
        FwEdgeIterator _e_out;

        ContractJob(const BwEdgeIterator& e_in, const NodeIterator& x, const FwEdgeIterator& e_out) noexcept
        : _e_in(e_in), _x(x), _e_out(e_out)
        {}
    };

    bool has_smaller_cost(const NodeIterator& u, const NodeIterator& v) const
    {
        assert( u < _nodes.size() );
        assert( v < _nodes.size() );
        assert( u < _random_permutation.size() );
        assert( v < _random_permutation.size() );

        if ( _node_cost[u] < _node_cost[v] )
            return true;

        if ( _node_cost[u] == _node_cost[v] )
            if ( _random_permutation[u] < _random_permutation[v] )
                return true;

        return false;
    }

    static bool dbg_check_merged_middle_node_data
    (
            const std::vector<MiddleNodeDescr> merged_middle_node_data,
            const EdgeData& present_edge,
            const TTF& new_ttf,
            const NodeIterator new_middle_node
    )
    {
        const size_t N_TEST_CASES = 20;

        for ( size_t i = 0 ; i != merged_middle_node_data.size() ; ++i )
        {
            const NodeIterator middle_node = merged_middle_node_data[i]._middle_node;

            const double a = merged_middle_node_data[i]._time;

            const double b =
                    i+1 < merged_middle_node_data.size() ?
                            merged_middle_node_data[i+1]._time :
                            merged_middle_node_data.front()._time + TTF::PERIOD();

            for ( size_t j = 0 ; j != N_TEST_CASES ; ++j )
            {
                const double x = util::random(a, b);

                if (
                    ! KATCH_EQUIV(
                            middle_node != new_middle_node,
                            middle_node == present_edge.get_middle_node(x)
                    )
                )
                    false;

                if (
                    ! KATCH_EQUIV(
                            middle_node == new_middle_node,
                            middle_node != present_edge.get_middle_node(x)
                    )
                )
                    false;

                if (
                    ! KATCH_IMPLIES(
                            present_edge.get_ttf().eval(x) + 0.0001 < new_ttf.eval(x),
                            middle_node != new_middle_node && present_edge.get_middle_node(x)
                    )
                )
                    false;

                if (
                    ! KATCH_IMPLIES(
                            present_edge.get_ttf().eval(x) > new_ttf.eval(x) + 0.0001,
                            middle_node == new_middle_node
                    )
                )
                    false;

                if (
                    ! KATCH_IMPLIES(
                            middle_node != new_middle_node,
                            present_edge.get_ttf().eval(x) <= new_ttf.eval(x) + 0.0001
                    )
                )
                    false;

                if (
                    ! KATCH_IMPLIES(
                            middle_node == new_middle_node,
                            present_edge.get_ttf().eval(x) + 0.0001 >= new_ttf.eval(x)
                    )
                )
                    false;
            }
        }

        return true;
    }

    void update_middle_node_data(EdgeData& edge, const TTF& new_ttf, const NodeIterator& new_middle_node)
    {
        std::deque<std::pair<double,RelativeTtfPosition>>
            relative_positioning = analyze_relative_position(new_ttf, edge.get_ttf());

        if ( relative_positioning.size() == 1 )
        {
            if ( relative_positioning.back().second == RelativeTtfPosition::G_ABOVE )
                edge.set_unique_middle_node(new_middle_node);

            return;
        }

        auto get_position =
                [&](const int& i)
                {
                    const int n(relative_positioning.size());
                    assert( ( i >= -1 && i <= n) );

                    double time;
                    RelativeTtfPosition position;

                    if ( i == -1 )
                    {
                        time = relative_positioning.back().first - TTF::PERIOD();
                        position = relative_positioning.back().second;
                    }
                    else if ( i == n )
                    {
                        time = relative_positioning.front().first + TTF::PERIOD();
                        position = relative_positioning.front().second;
                    }
                    else
                    {
                        time = relative_positioning[i].first;
                        position = relative_positioning[i].second;
                    }

                    return std::make_pair(time, position);
                };

        size_t i = 0;
        size_t j = 0;

        std::vector<MiddleNodeDescr> middle_node_data;

        while ( i < relative_positioning.size() || j < edge.get_n_middle_nodes() )
        {
            if ( eq( get_position(i).first, edge.get_middle_node_descr(j)._time ) )
            {
                if ( get_position(i).second == RelativeTtfPosition::G_ABOVE )
                    middle_node_data.emplace_back( get_position(i).first, new_middle_node );
                else
                    middle_node_data.push_back( edge.get_middle_node_descr(j) );

                ++i; ++j;
            }
            else if ( get_position(i).first < edge.get_middle_node_descr(j)._time )
            {
                if ( get_position(i).second == RelativeTtfPosition::G_ABOVE )
                    middle_node_data.emplace_back( get_position(i).first, new_middle_node );
                else
                    middle_node_data.push_back(
                            MiddleNodeDescr( get_position(i).first, edge.get_middle_node_descr(j-1)._middle_node) );

                ++i;
            }
            else
            {
                assert( edge.get_middle_node_descr(j)._time < get_position(i).first );

                if ( get_position(i-1).second == RelativeTtfPosition::F_ABOVE )
                    middle_node_data.push_back( edge.get_middle_node_descr(j) );

                ++j;
            }
        }

        assert( dbg_check_merged_middle_node_data(middle_node_data, edge, new_ttf, new_middle_node) );

        edge.set_middle_node_data(middle_node_data.begin(), middle_node_data.end());
    }

    bool is_minimum_in_neighborhood(const NodeIterator& x, const int& hop_radius) const
    {
        assert( x < _graph.get_n_nodes() );
        assert( ! has_smaller_cost(x, x) );

        std::queue<NodeIterator> Q;
        Q.push(x);

        std::unordered_map<NodeIterator,int> n_hops;
        n_hops[x] = 0;

        while ( ! Q.empty() )
        {
            const NodeIterator u = Q.front();
            Q.pop();

            for ( FwEdgeIterator e = _graph.out_edges_begin(u) ; e != _graph.out_edges_end(u) ; ++e )
            {
                const NodeIterator v = _graph.get_target(e);
                if ( n_hops.find(v) != n_hops.end() ) continue;

                assert( has_smaller_cost(v, x) || has_smaller_cost(x, v) );
                if ( has_smaller_cost(v, x) ) return false;

                n_hops[v] = n_hops[u] + 1;

                if ( n_hops[v] >= hop_radius ) continue;
                Q.push(v);
            }

            for ( BwEdgeIterator e = _graph.in_edges_begin(u) ; e != _graph.in_edges_end(u) ; ++e )
            {
                const NodeIterator v = _graph.get_source(e);
                if ( n_hops.find(v) != n_hops.end() ) continue;

                assert( has_smaller_cost(v, x) || has_smaller_cost(x, v) );
                if ( has_smaller_cost(v, x) ) return false;

                n_hops[v] = n_hops[u] + 1;

                if ( n_hops[v] >= hop_radius ) continue;
                Q.push(v);
            }
}

        return true;
    }

    // Choose nodes to be contracted next and move them in <code>_nodes</code>
    // such that they occur after <code>_working_nodes_begin</code>. This includes
    // setting up <code>_working_nodes_end</code>, such that
    //
    //   [_working_nodes_begin, _working_nodes_end)
    //
    // contains the nodes to be contracted next afterwards.
    void parallelly_choose_nodes_to_be_contracted_next()
    {
        assert( _working_nodes_end <= _nodes.size() );
        assert( ! util::has_duplicates(&_nodes[_working_nodes_end], &_nodes[_nodes.size()]) );

        // chose nodes to be contracted next
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
            for ( size_t i = _working_nodes_end ; i < _nodes.size() ; ++i )
            {

                const NodeIterator x = _nodes[i];
                assert( ! _node_selected_for_contraction[x] );
                const bool x_is_local_minimum = is_minimum_in_neighborhood(x, 2);

                if ( x_is_local_minimum )
                {
                    assert( x < _node_selected_for_contraction.size() );
                    assert( _node_selected_for_contraction[x] == 0 );
                    _node_selected_for_contraction[x] = 1;
                }
            }
       }

        // move these nodes to the right position in <code>_nodes</code>
        auto first_unselected_node_it =
                std::stable_partition(_nodes.begin() + _working_nodes_end, _nodes.end(),
                       [&](const NodeIterator& x) { return _node_selected_for_contraction[x]; } );

        _working_nodes_end = first_unselected_node_it - _nodes.begin();

        assert( _working_nodes_begin < _working_nodes_end );
        assert( _working_nodes_end <= _nodes.size() );
        assert( ! util::has_duplicates(&_nodes[_working_nodes_begin], &_nodes[_working_nodes_end]) );
        assert( ! util::has_duplicates(_nodes.begin(), _nodes.end()) );
        assert( ! util::has_duplicates( util::neighborhood(_graph, &_nodes[_working_nodes_begin], &_nodes[_working_nodes_end]) ) );
    }


    bool is_shortcut_necessary
    (
            const BwEdgeIterator& e_in,
            const NodeIterator& x,
            const FwEdgeIterator& e_out,
            const bool simulate,
            EdgeDescription& insert_edge
    ) const
    {
        const NodeIterator u = _graph.get_source(e_in);
        const NodeIterator v = _graph.get_target(e_out);
        const EdgeData& edge_ux = _graph.get_edge_data(e_in);
        const EdgeData& edge_xv = _graph.get_edge_data(e_out);

        assert( u != v );
        assert( u != x );
        assert( v != x );

        // use already cached witness if present
        WitnessCacheEntry witness;
        assert( witness._shortcut_status == WitnessSearch::UNDECIDED );

        if ( _witness_cache.contains(u,x,v) )
        {
            witness = _witness_cache.lookup(u,x,v);
            assert( witness._u == u  && witness._x == x && witness._v == v );
            assert( witness._shortcut_status != WitnessSearch::UNDECIDED );
        }

        size_t n_original_edges = edge_ux.get_n_original_edges() + edge_xv.get_n_original_edges();

        size_t shortcut_complexity =
                    witness._shortcut_status == WitnessSearch::UNDECIDED ?
                                  (edge_ux.get_ttf().size() + edge_xv.get_ttf().size()) :
                                  witness._shortcut_complexity;

        insert_edge = EdgeDescription(u, v, x, shortcut_complexity, n_original_edges);

        // return if witness is already cached
        if ( witness._shortcut_status != WitnessSearch::UNDECIDED )
                return witness._shortcut_status == WitnessSearch::NECESSARY;

        assert( omp_get_thread_num() < int( _local_thread_data.size() ) );
        LocalThreadData& thread_data = _local_thread_data[omp_get_thread_num()];

        // perform witness search if witness is not cached yet
        const ShortcutStatus shortcut_status = thread_data._witness_search.run(u, edge_ux, x, edge_xv, v, insert_edge);
        assert( shortcut_status != WitnessSearch::UNDECIDED );

        if ( simulate )
            thread_data._witnesses_to_insert_to_cache.emplace_back(shortcut_status, insert_edge._complexity, u, x, v);

        return shortcut_status == WitnessSearch::NECESSARY;
    }

    void simulate_contraction_and_update_cost_of_node(NodeIterator x)
    {
        uint32_t n_insert_edges = 0;

        uint32_t n_original_edges_of_inserted_edges = 0;
        uint32_t n_original_edges_of_removed_edges = 0;

        uint32_t total_complexity_of_inserted_edges = 0;
        uint32_t total_complexity_of_removed_edges = 0;

        uint32_t n_removed_edges = 0;

        for ( FwEdgeIterator e = _graph.out_edges_begin(x) ; e != _graph.out_edges_end(x) ; ++e )
        {
            ++n_removed_edges;
            n_original_edges_of_removed_edges += _graph.get_edge_data(e).get_n_original_edges();
            total_complexity_of_removed_edges += _graph.get_edge_data(e).get_ttf().size();
        }

        for ( BwEdgeIterator e = _graph.in_edges_begin(x) ; e != _graph.in_edges_end(x) ; ++e )
        {
            ++n_removed_edges;
            n_original_edges_of_removed_edges += _graph.get_edge_data(e).get_n_original_edges();
            total_complexity_of_removed_edges += _graph.get_edge_data(e).get_ttf().size();
        }

        for ( BwEdgeIterator e_in = _graph.in_edges_begin(x) ; e_in != _graph.in_edges_end(x) ; ++e_in )
        {
            for ( FwEdgeIterator e_out = _graph.out_edges_begin(x) ; e_out != _graph.out_edges_end(x) ; ++e_out )
            {
                if ( _graph.get_source(e_in) == _graph.get_target(e_out) )
                    continue;

                EdgeDescription insert_edge;
                bool shortcut_necessary = is_shortcut_necessary(e_in, x, e_out, SIMULATE, insert_edge);

                if ( shortcut_necessary )
                {
                    ++n_insert_edges;

                    n_original_edges_of_inserted_edges += insert_edge._n_original_edges;
                    total_complexity_of_inserted_edges += insert_edge._complexity;
                }
            }
        }

        double edges_quotient = double(n_insert_edges) / std::max(1.0, double(n_removed_edges));

        double original_edges_quotient =
                double(n_original_edges_of_inserted_edges) / std::max(1.0, double(n_original_edges_of_removed_edges));

        double complexity_quotient =
                double(total_complexity_of_inserted_edges) / std::max(1.0, double(total_complexity_of_removed_edges));

        double new_node_cost =
                    2.0 * edges_quotient +
                    _hierarchie_depth[x] +
                    original_edges_quotient +
                    2.0 * complexity_quotient;

        _node_cost[x] = new_node_cost;
    }
    
    void contract(const BwEdgeIterator& e_in, const NodeIterator x, const FwEdgeIterator& e_out, bool simulate = DO_NOT_SIMULATE)
    {
        EdgeDescription insert_edge;
        const bool shortcut_necessary = is_shortcut_necessary(e_in, x, e_out, simulate, insert_edge);

        if ( ! shortcut_necessary ) return;

        if ( insert_edge._ttf == TTF::INFTY )
        {
            const EdgeData& edge_ux = _graph.get_edge_data(e_in);
            const EdgeData& edge_xv = _graph.get_edge_data(e_out);

            insert_edge._ttf = link(edge_xv.get_ttf(), edge_ux.get_ttf());
        }

        assert( omp_get_thread_num() < int(_local_thread_data.size()) );
        _local_thread_data[omp_get_thread_num()]._edges_to_insert.push_back(std::move(insert_edge));

    }

    void parallelly_contract_nodes()
    {
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
            for ( size_t i = _working_nodes_begin ; i < _working_nodes_end ; ++i )
            {
                NodeIterator x = _nodes[i];

                for ( BwEdgeIterator e_in = _graph.in_edges_begin(x) ; e_in != _graph.in_edges_end(x) ; ++e_in )
                    for ( FwEdgeIterator e_out = _graph.out_edges_begin(x) ; e_out != _graph.out_edges_end(x) ; ++e_out )
                        if ( _graph.get_source(e_in) != _graph.get_target(e_out) )
                            contract(e_in, x, e_out);
            }
        }

        for ( size_t i = _working_nodes_begin ; i < _working_nodes_end ; ++i )
        {
            NodeIterator x = _nodes[i];
            _witness_cache.remove(x);

            for ( BwEdgeIterator e_in = _graph.in_edges_begin(x) ; e_in != _graph.in_edges_end(x) ; ++e_in )
            {
                NodeIterator u = _graph.get_source(e_in);

                for ( FwEdgeIterator e_out = _graph.out_edges_begin(x) ; e_out != _graph.out_edges_end(x) ; ++e_out )
                {
                    NodeIterator v = _graph.get_target(e_out);

                    _witness_cache.remove(INVALID_NODE_ITERATOR, u, x);
                    _witness_cache.remove(x, v, INVALID_NODE_ITERATOR);
                }
            }
        }

        #ifndef NDEBUG
            for ( auto& thread_data : _local_thread_data )
                assert( thread_data._witnesses_to_insert_to_cache.empty() );
        #endif

    }


    // the graph to be transformed into a TCH
    Graph _graph;

    // [0, _working_nodes_begin)                  : already contracted nodes
    // [_working_nodes_begin, _working_nodes_end) : nodes to be contracted next
    // [_working_nodes_end, nodes.size())         : not yet contracted nodes
    size_t _working_nodes_begin;
    size_t _working_nodes_end;
    std::vector<NodeIterator> _nodes;

    std::vector<char> _node_selected_for_contraction;

    // maps nodes to their tentative cost
    std::vector<double> _node_cost;

    // random permutation tie breaking
    std::vector<int> _random_permutation;

    // part of the tentative node cost
    std::vector<int> _hierarchie_depth;

    // maps node iterator to level
    std::vector<uint32_t> _node_to_level;

    // cache for the results of witness searches
    mutable WitnessCache _witness_cache;

    // file where the tch is stored to
    btch_format::OutputFile _btch_output_file;

    // local data of the parallel threads
    mutable std::vector<LocalThreadData> _local_thread_data;

    // the value added to sample search score when a shortcut is added
    const int SAMPLE_SEARCH_LAMBDA_PLUS = 4;

    // the value subtracted from sample search score when a shortcut is omitted
    const int SAMPLE_SEARCH_LAMBDA_MINUS = 1;

    // the threshold for turning sample search mode on (if score exceeds XI) and off (if score gets smaller than -XI)
    const int SAMPLE_SEARCH_XI = 1000;


public:

    Ordering(Graph&& graph, const std::string& tch_outout_file_name)
    : _graph(std::move(graph)),
      _working_nodes_begin(0), _working_nodes_end(0), _nodes(_graph.get_n_nodes(), INVALID_NODE_ITERATOR),
      _node_selected_for_contraction(_graph.get_n_nodes(), 0),
      _node_cost(_graph.get_n_nodes(), 0.0),
      _random_permutation(_graph.get_n_nodes(), INVALID_NODE_ITERATOR),
      _hierarchie_depth(_graph.get_n_nodes(), 1),
      _node_to_level(_graph.get_n_nodes(), INVALID_NODE_ITERATOR),
      _witness_cache(_graph.get_n_nodes()),
      _btch_output_file(tch_outout_file_name, _graph.get_n_nodes(), TTF::PERIOD()),
      _local_thread_data()
    {
        // set up node set
        std::generate(_nodes.begin(), _nodes.end(), []() { static size_t counter(0); return counter++; } );

        // set up random permutation used for tie breaking when tentative node costs are compared
        std::generate(_random_permutation.begin(), _random_permutation.end(), []() { static size_t counter(0); return counter++; } );
        std::random_shuffle(_random_permutation.begin(), _random_permutation.end());

        assert( ! util::has_duplicates(_nodes) );
        assert( ! util::has_duplicates(_random_permutation) );
    }

    void order_and_construct(int n_threads = MAXIMUM_THREADS)
    {
        auto t_begin = util::time_stamp();

        double total_io_time = 0;
        _local_thread_data.clear();

        if ( n_threads == MAXIMUM_THREADS )
            omp_set_num_threads( omp_get_max_threads() );
        else
            omp_set_num_threads( std::min(n_threads, omp_get_max_threads()) );

        _local_thread_data.reserve(omp_get_max_threads());
        for ( int i = 0 ; i < omp_get_max_threads() ; ++i )
            _local_thread_data.emplace_back(&_graph);

        KATCH_STATUS("Preprocessing running with " << std::min(n_threads, omp_get_max_threads()) << " threads\n");

        //
        // set up initial tentative node costs
        //
        KATCH_STATUS("Computing initial node costs... ");
        auto t1 = util::time_stamp();
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
            for ( NodeIterator x = _graph.nodes_begin() ; x < _graph.nodes_end() ; ++x )
                simulate_contraction_and_update_cost_of_node(x);
        }

        //
        // update witness cache
        //
        for ( auto& thread_data : _local_thread_data )
        {
            for ( auto witness_cache_entry : thread_data._witnesses_to_insert_to_cache )
                _witness_cache.insert(witness_cache_entry);

            thread_data._witnesses_to_insert_to_cache.clear();
        }

        auto t2 = util::time_stamp();
        KATCH_CONTINUE_STATUS(katch::util::get_duration_in_seconds(t1, t2) << " sec\n");

        //
        // repeatedly compute independent sets and contract their nodes,
        // merging necessary shortcuts with already present edges
        //
        KATCH_STATUS("Contract independent node sets of size ");
        while ( _working_nodes_begin < _nodes.size() )
        {
            //
            // move the nodes to be contracted next to appear after <code>_working_nodes_begin</code>
            // in <code>_nodes</code> and set up <code>_working_nodes_end</code>
            //
            parallelly_choose_nodes_to_be_contracted_next();

            assert( _working_nodes_begin < _working_nodes_end );
            assert( _working_nodes_end <= _nodes.size() );
            assert( ! util::has_duplicates(&_nodes[_working_nodes_begin], &_nodes[_working_nodes_end]) );
            KATCH_CONTINUE_STATUS( _working_nodes_end - _working_nodes_begin );

            //
            // perform witness search for the nodes to contracted and store the new and to be merged shortcuts
            //
            parallelly_contract_nodes();

            KATCH_CONTINUE_STATUS(" ");

            //
            // store incident edges of the contracted nodes to the .btch-file
            // also remember the nodes adjacent to the contracted nodes to update their cost later
            // also remove all incident edges of the contracted nodes
            //
            std::vector<NodeIterator> nodes_to_update;
            for ( size_t i = _working_nodes_begin ; i < _working_nodes_end ; ++i )
            {
                assert( i < _nodes.size() );
                const NodeIterator x = _nodes[i];

                assert( x < _node_to_level.size() );
                assert( _node_to_level[x] == INVALID_NODE_ITERATOR );
                _node_to_level[x] = i;

                for ( FwEdgeIterator e = _graph.out_edges_begin(x) ; e != _graph.out_edges_end(x) ; ++e )
                {
                    nodes_to_update.push_back( _graph.get_target(e) );

                    _hierarchie_depth[_graph.get_target(e)] =
                            std::max(_hierarchie_depth[_graph.get_target(e)], _hierarchie_depth[x] + 1);

                    auto t3 = util::time_stamp();
                    _btch_output_file.write_edge(x, _graph.get_target(e), _graph.get_edge_data(e));
                    auto t4 = util::time_stamp();
                    total_io_time += util::get_duration_in_seconds(t3, t4);
                }

                for ( BwEdgeIterator e = _graph.in_edges_begin(x) ; e != _graph.in_edges_end(x) ; ++e )
                {
                    nodes_to_update.push_back( _graph.get_source(e) );

                    _hierarchie_depth[_graph.get_source(e)] =
                            std::max(_hierarchie_depth[_graph.get_source(e)], _hierarchie_depth[x] + 1);

                    auto t3 = util::time_stamp();
                    _btch_output_file.write_edge(_graph.get_source(e), x, _graph.get_edge_data(e));
                    auto t4 = util::time_stamp();
                    total_io_time += util::get_duration_in_seconds(t3, t4);
                }

                _graph.remove_all_edges_of_node(x);
            }

            assert( _graph.dbg_check_graph() );

            util::remove_duplicates(nodes_to_update);

            KATCH_CONTINUE_STATUS("-");

            //
            // insert/merge new shortcut edges and remove merged edges from cache (otherwise computation will be wrong)
            //
            for ( auto& thread_data : _local_thread_data )
            {
                for ( auto& edge_descr : thread_data._edges_to_insert )
                {
                    assert( edge_descr._ttf != TTF::INFTY );
                    FwEdgeIterator e = std::get<0>(_graph.get_edge_iterators(edge_descr._src, edge_descr._tgt));

                    if ( e == INVALID_FW_EDGE_ITERATOR )
                    {
                        EdgeData edge_data;
                        edge_data.set_ttf( std::move(edge_descr._ttf) );
                        edge_data.set_n_original_edges(edge_descr._n_original_edges);
                        edge_data.set_unique_middle_node(edge_descr._middle_node);

                        _graph.insert_edge(edge_descr._src, edge_descr._tgt, std::move(edge_data));
                    }
                    else
                    {
                        EdgeData& edge_data = _graph.get_edge_data(e);

                        update_middle_node_data(edge_data, edge_descr._ttf, edge_descr._middle_node);
                        edge_data.set_ttf( merge(edge_data.get_ttf(), edge_descr._ttf) );

                        _witness_cache.remove(INVALID_NODE_ITERATOR, edge_descr._src, edge_descr._tgt);
                        _witness_cache.remove(edge_descr._src, edge_descr._tgt, INVALID_NODE_ITERATOR);
                    }
                }

                thread_data._edges_to_insert.clear();
            }

            assert( _graph.dbg_check_graph() );

            KATCH_CONTINUE_STATUS(" ");

            //
            // update cost of nodes adjacent to contracted nodes
            //
            #pragma omp parallel
            {
                #pragma omp for schedule(dynamic)
                for ( size_t i = 0 ; i < nodes_to_update.size() ; ++i )
                {
                    NodeIterator x = nodes_to_update[i];
                    simulate_contraction_and_update_cost_of_node(x);
                }
            }

            //
            // update witness cache
            //
            for ( auto& thread_data : _local_thread_data )
            {
                for ( auto witness_cache_entry : thread_data._witnesses_to_insert_to_cache )
                    _witness_cache.insert(witness_cache_entry);

                thread_data._witnesses_to_insert_to_cache.clear();
            }


            //
            // mark all nodes in the current independent set as contracted
            //
            _working_nodes_begin = _working_nodes_end;

        }

        KATCH_CONTINUE_STATUS("finished\n");

        assert( _node_to_level.size() == _nodes.size() );

        assert( std::all_of(_node_to_level.begin(), _node_to_level.end(),
                [&](const uint32_t& level) -> bool { return size_t(level) < _graph.get_n_nodes(); }) );

        assert( ! util::has_duplicates(_node_to_level) );


        auto t7 = util::time_stamp();
        KATCH_STATUS("Write node level information to BTCH file... ");
        _btch_output_file.write_level_info(_node_to_level);
        KATCH_CONTINUE_STATUS("OK\n");
        KATCH_STATUS("Finalize BTCH file... ");
        _btch_output_file.close();
        KATCH_CONTINUE_STATUS("OK\n");
        auto t8 = util::time_stamp();
        total_io_time += util::get_duration_in_seconds(t7, t8);

        auto t_end = util::time_stamp();
        double total_ordering_time = katch::util::get_duration_in_seconds(t_begin, t_end);

        KATCH_STATUS("Ordering took " << total_ordering_time << " sec in total.\n");
        KATCH_STATUS("Writing hierarchy to file took " << total_io_time << " sec in total.\n");
        KATCH_STATUS("Ordering took " << (total_ordering_time - total_io_time) << " sec without writing to file.\n");
    }
};

const Ordering::WitnessCacheEntry Ordering::WitnessCache::INVALID_ENTRY {};

}

#endif /* KATCH_PREPROCESSING_H_ */
