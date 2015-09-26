/*
 * katch/datastr/graph/dynamic_graph.h
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
 * License along with KaTCH; see the file COPYING; if not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef KATCH_DYNAMIC_GRAPH_H_
#define KATCH_DYNAMIC_GRAPH_H_

#include <algorithm>
#include <assert.h>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>
#include <limits>

#include "datastr/graph/basic.h"
#include "datastr/graph/order_edge_data.h"
#include "util/misc.h"


namespace katch
{

template <typename EdgeDataT>
class DynamicGraph
{

private:

    static constexpr double GROWTH_FACTOR = 1.2;

public:

    using EdgeData = EdgeDataT;

    struct ShortcutDataItem
    {
        double _time;
        NodeIterator _middle_node;
    };

private:

    class Node
    {
    public:

        Node()
        : _out_edges_begin(INVALID_FW_EDGE_ITERATOR),
          _out_edges_end(INVALID_BW_EDGE_ITERATOR),
          _in_edges_begin(INVALID_FW_EDGE_ITERATOR),
          _in_edges_end(INVALID_FW_EDGE_ITERATOR)
        {}

        FwEdgeIterator _out_edges_begin;
        FwEdgeIterator _out_edges_end;
        BwEdgeIterator _in_edges_begin;
        BwEdgeIterator _in_edges_end;
    };

    struct ForwardEdge
    {
        NodeIterator _target;
        BwEdgeIterator _backward_edge;
        EdgeData _data;

        ForwardEdge(const ForwardEdge&) = delete;
        ForwardEdge& operator= (const ForwardEdge&) = delete;

        ForwardEdge(ForwardEdge&&) = default;
        ForwardEdge& operator= (ForwardEdge&&) = default;

        ForwardEdge()
        : _target(INVALID_NODE_ITERATOR), _backward_edge(INVALID_BW_EDGE_ITERATOR), _data()
        {
            assert( is_dummy() );
        }

        template <typename EdgeInfo>
        void set_up(EdgeInfo&& edge_info)
        {
            _target = edge_info.get_target();
            _data.set_up(std::move(edge_info));
        }

        // tells whether this edge is a dummy entry
        bool is_dummy() const
        {
            return _target == INVALID_NODE_ITERATOR;
        }

        // make this edge a dummy entry
        void make_dummy()
        {
            _target = INVALID_NODE_ITERATOR;
            _data = std::move(EdgeData());
        }
    };

    struct BackwardEdge
    {
        NodeIterator _source;
        FwEdgeIterator _forward_edge;

        BackwardEdge()
        : _source(INVALID_NODE_ITERATOR), _forward_edge(INVALID_FW_EDGE_ITERATOR)
        {
            assert( is_dummy() );
        }

        // tells whether this edge is a dummy entry
        bool is_dummy() const
        {
            return _source == INVALID_NODE_ITERATOR;
        }

        // make this edge a dummy entry
        void make_dummy()
        {
            _source = INVALID_NODE_ITERATOR;
        }
    };

    void remove_out_edge(const NodeIterator& src, const FwEdgeIterator& e)
    {
        assert( _nodes[src]._out_edges_begin <= e );
        assert( e < _nodes[src]._out_edges_end );
        assert( e < FwEdgeIterator( _forward_edges.size()) );
        assert( ! _forward_edges[e].is_dummy() );

        assert( _nodes[src]._out_edges_end > 0 );
        FwEdgeIterator last_edge_it( _nodes[src]._out_edges_end - 1 );
        assert( last_edge_it < _forward_edges.size() );

        if ( e != last_edge_it )
        {
            _forward_edges[e] = std::move( _forward_edges[last_edge_it] );
            _backward_edges[_forward_edges[e]._backward_edge]._forward_edge = e;
        }

        _forward_edges[last_edge_it].make_dummy();
        --_nodes[src]._out_edges_end;
    }

    void remove_in_edge(const NodeIterator& tgt, const BwEdgeIterator& e)
    {
        assert( _nodes[tgt]._in_edges_begin <= e );
        assert( e < _nodes[tgt]._in_edges_end );
        assert( e < BwEdgeIterator( _backward_edges.size()) );
        assert( ! _backward_edges[e].is_dummy() );

        assert( _nodes[tgt]._in_edges_end > 0 );
        BwEdgeIterator last_edge_it( _nodes[tgt]._in_edges_end - 1 );
        assert( last_edge_it < _backward_edges.size() );
        assert( ! _backward_edges[last_edge_it].is_dummy() );

        if ( e != last_edge_it )
        {
            _backward_edges[e] = _backward_edges[last_edge_it];
            _forward_edges[_backward_edges[e]._forward_edge]._backward_edge = e;
        }

        _backward_edges[last_edge_it].make_dummy();
        --_nodes[tgt]._in_edges_end;
    }

    static inline size_t compute_n_dummies(const size_t& k)
    {
        size_t result = std::max( k + 2, size_t(double(k) * GROWTH_FACTOR) ) - k;

        assert( result >= 2 );
        return result;
    }

    inline BwEdgeIterator& access_iterator_entry_in_reverse_edge(const BwEdgeIterator& e)
    {
        assert( e == _forward_edges[_backward_edges[e]._forward_edge]._backward_edge );
        return _forward_edges[_backward_edges[e]._forward_edge]._backward_edge;
    }

    inline FwEdgeIterator& access_iterator_entry_in_reverse_edge(const FwEdgeIterator& e)
    {
        assert( e == _backward_edges[_forward_edges[e]._backward_edge]._forward_edge );
        return _backward_edges[_forward_edges[e]._backward_edge]._forward_edge;
    }

    inline static bool is_invalid_edge_iterator(const BwEdgeIterator& e)
    {
        return e == INVALID_BW_EDGE_ITERATOR;
    }

    inline static bool is_invalid_edge_iterator(FwEdgeIterator e)
    {
        return e == INVALID_FW_EDGE_ITERATOR;
    }

    template <typename EdgeType, typename EdgeIteratorType>
    EdgeIteratorType insert_edge_impl(std::vector<EdgeType>& edges, EdgeIteratorType& begin, EdgeIteratorType& end)
    {
        EdgeIteratorType edge_insert_it;

        assert( KATCH_EQUIV( is_invalid_edge_iterator(begin),  is_invalid_edge_iterator(end) ) );

        if ( is_invalid_edge_iterator(end) )
        {
            const size_t n_dummies = compute_n_dummies(1) + 1;
            for ( size_t i = 0 ; i < n_dummies ; ++i ) edges.emplace_back();

            begin = edges.size();
            edge_insert_it = edges.size();
            end = edges.size() + 1;
        }
        else if ( end == edges.size() )
        {
            const size_t old_n_edges = end - begin;
            const size_t n_dummies = compute_n_dummies(1 + old_n_edges) + 1;
            for ( size_t i = 0 ; i < n_dummies ; ++i ) edges.emplace_back();

            edge_insert_it = end;
            ++end;
        }
        else if ( edges[end].is_dummy() )
        {
            edge_insert_it = end;
            ++end;
        }
        else if ( begin > 0 && edges[begin - 1].is_dummy() )
        {
            --begin;
            edge_insert_it = begin;
        }
        else
        {
            const size_t old_n_edges = end - begin;
            EdgeIteratorType begin_new(edges.size());

            for ( EdgeIteratorType it = begin ; it != end ; ++it )
            {
                access_iterator_entry_in_reverse_edge(it) = edges.size();
                edges.emplace_back( std::move(edges[it]) );
                edges[it].make_dummy();
            }

            begin = begin_new;
            edge_insert_it = edges.size();
            end = edges.size() + 1;

            const size_t n_dummies = compute_n_dummies(1 + old_n_edges) + 1;
            for ( size_t i = 0 ; i < n_dummies ; ++i ) edges.emplace_back();
        }

        assert( ! is_invalid_edge_iterator(edge_insert_it) );
        assert( edge_insert_it < edges.size() );
        assert( edges[edge_insert_it].is_dummy() );
        return edge_insert_it;
    }

    std::vector<Node> _nodes;
    std::vector<ForwardEdge> _forward_edges;
    std::vector<BackwardEdge> _backward_edges;
    size_t _n_edges;

public:

    DynamicGraph()
    : _nodes(), _forward_edges(), _backward_edges(), _n_edges(0)
    {}

    DynamicGraph(const DynamicGraph&) = delete;
    DynamicGraph& operator= (const DynamicGraph&) = delete;

    DynamicGraph(DynamicGraph&& graph) = default;
    DynamicGraph& operator= (DynamicGraph&& graph) = default;

    template <typename EdgeInfo>
    DynamicGraph(std::vector<EdgeInfo>&& edge_list)
    : _nodes(), _forward_edges(), _backward_edges(), _n_edges(edge_list.size())
    {
        size_t max_node_it = 0;

        std::for_each( edge_list.begin(), edge_list.end(),
                        [&max_node_it](const EdgeInfo& edge_info) -> void
                        {
                            max_node_it = std::max(max_node_it, edge_info.get_source());
                            max_node_it = std::max(max_node_it, edge_info.get_target());
                        } );

        _nodes.assign(max_node_it + 1, Node());

        //
        // set up forward edges
        //

        {
            std::sort( edge_list.begin(), edge_list.end(),
                            [](const EdgeInfo& lhs, const EdgeInfo& rhs) -> bool
                                    { return lhs.get_source() < rhs.get_source(); } );

            for ( auto it = edge_list.begin() ; it != edge_list.end() ; ++it )
            {
                if ( _nodes[it->get_source()]._out_edges_begin == INVALID_FW_EDGE_ITERATOR )
                {
                    // insert dummy edges
                    if ( it > edge_list.begin() )
                    {
                        NodeIterator previous_source = (it-1)->get_source();
                        size_t n_edges_of_previos_source = _nodes[previous_source]._out_edges_end - _nodes[previous_source]._out_edges_begin;
                        size_t n_dummies = std::max( size_t(n_edges_of_previos_source + 2), size_t(double(n_edges_of_previos_source) * GROWTH_FACTOR) );

                        for ( size_t i = 0 ; i < n_dummies ; ++i )
                            _forward_edges.push_back(std::move(ForwardEdge()));
                    }

                    _nodes[it->get_source()]._out_edges_begin = _forward_edges.size();
                    _nodes[it->get_source()]._out_edges_end = _forward_edges.size();
                }

                _forward_edges.push_back(ForwardEdge());
                _forward_edges.back().set_up( std::move(*it) );

                ++_nodes[it->get_source()]._out_edges_end;
            }

            if ( ! edge_list.empty() )
            {
                NodeIterator previous_source = edge_list.back().get_source();
                size_t n_edges_of_previos_source = _nodes[previous_source]._out_edges_end - _nodes[previous_source]._out_edges_begin;
                size_t n_dummies = std::max( size_t(n_edges_of_previos_source + 2), size_t(double(n_edges_of_previos_source) * GROWTH_FACTOR) );

                for ( size_t i = 0 ; i < n_dummies ; ++i )
                    _forward_edges.push_back(std::move(ForwardEdge()));
            }

            edge_list.clear();
        }

        //
        // set up backward edges
        //

        {
            struct BackwardEdgeInfo
            {
                NodeIterator _source;
                NodeIterator _target;
                FwEdgeIterator _forward_edge;
            };

            std::vector<BackwardEdgeInfo> backward_edge_list;
            backward_edge_list.reserve(_n_edges);

            for ( NodeIterator source = 0 ; source != _nodes.size() ; ++source )
            {
                for ( FwEdgeIterator fw_edge_it = _nodes[source]._out_edges_begin ; fw_edge_it < _nodes[source]._out_edges_end ; ++fw_edge_it )
                {
                    assert( fw_edge_it < _forward_edges.size() );

                    backward_edge_list.push_back(BackwardEdgeInfo());
                    backward_edge_list.back()._source = source;
                    backward_edge_list.back()._target = _forward_edges[fw_edge_it]._target;
                    backward_edge_list.back()._forward_edge = fw_edge_it;
                }
            }

            std::sort(
                backward_edge_list.begin(), backward_edge_list.end(),
                    [](const BackwardEdgeInfo& a, const BackwardEdgeInfo& b) -> bool { return a._target < b._target; }
            );

            for ( auto it = backward_edge_list.begin() ; it != backward_edge_list.end() ; ++it )
            {
                assert( it->_target < _nodes.size() );

                if ( _nodes[it->_target]._in_edges_begin == INVALID_BW_EDGE_ITERATOR )
                {
                    // insert dummy edges
                    if ( it > backward_edge_list.begin() )
                    {
                        NodeIterator previous_target = (it-1)->_target;
                        size_t n_edges_of_previos_target = _nodes[previous_target]._in_edges_end - _nodes[previous_target]._in_edges_begin;
                        size_t n_dummies = std::max( size_t(n_edges_of_previos_target + 2), size_t(double(n_edges_of_previos_target) * GROWTH_FACTOR) );

                        _backward_edges.insert(_backward_edges.end(), n_dummies, BackwardEdge());
                    }

                    _nodes[it->_target]._in_edges_begin = _backward_edges.size();
                    _nodes[it->_target]._in_edges_end = _backward_edges.size();
                }

                assert( it->_forward_edge < _forward_edges.size() );
                _forward_edges[it->_forward_edge]._backward_edge = _backward_edges.size();

                _backward_edges.push_back(BackwardEdge());
                _backward_edges.back()._source = it->_source;
                _backward_edges.back()._forward_edge = it->_forward_edge;

                ++_nodes[it->_target]._in_edges_end;
            }

            if ( ! backward_edge_list.empty() )
            {
                NodeIterator previous_target = backward_edge_list.back()._target;
                size_t n_edges_of_previos_source = _nodes[previous_target]._in_edges_end - _nodes[previous_target]._in_edges_begin;
                size_t n_dummies = std::max( size_t(n_edges_of_previos_source + 2), size_t(double(n_edges_of_previos_source) * GROWTH_FACTOR) );

                for ( size_t i = 0 ; i < n_dummies ; ++i )
                    _backward_edges.push_back(BackwardEdge());
            }

            backward_edge_list.clear();
        }

        assert( dbg_check_graph() );
        assert( ! has_parallel_edges() );
    }

    bool dbg_check_graph() const
    {
        for ( NodeIterator u = 0 ; u != _nodes.size() ; ++u )
        {
            assert(
                KATCH_EQUIV(
                    _nodes[u]._out_edges_begin == INVALID_FW_EDGE_ITERATOR,
                    _nodes[u]._out_edges_end == INVALID_FW_EDGE_ITERATOR
                )
            );

            assert( _nodes[u]._out_edges_begin == INVALID_FW_EDGE_ITERATOR ||
                                _nodes[u]._out_edges_begin <= _forward_edges.size() );
            assert( _nodes[u]._out_edges_end == INVALID_FW_EDGE_ITERATOR ||
                                _nodes[u]._out_edges_end <= _forward_edges.size() );

            for ( FwEdgeIterator e = _nodes[u]._out_edges_begin ; e != _nodes[u]._out_edges_end ; ++e )
            {
                assert( e != INVALID_FW_EDGE_ITERATOR );
                assert( e < _forward_edges.size() );
                assert( ! _forward_edges[e].is_dummy() );
                assert( _forward_edges[e]._backward_edge != INVALID_BW_EDGE_ITERATOR );
                assert( ! _backward_edges[_forward_edges[e]._backward_edge].is_dummy() );
                assert( _backward_edges[_forward_edges[e]._backward_edge]._forward_edge == e );
            }

            assert(
                KATCH_EQUIV(
                    _nodes[u]._in_edges_begin == INVALID_BW_EDGE_ITERATOR,
                    _nodes[u]._in_edges_begin == INVALID_BW_EDGE_ITERATOR
                )
            );
            assert( _nodes[u]._in_edges_begin == INVALID_BW_EDGE_ITERATOR || _nodes[u]._in_edges_begin < _backward_edges.size() );
            assert( _nodes[u]._in_edges_end == INVALID_BW_EDGE_ITERATOR || _nodes[u]._in_edges_end < _backward_edges.size() );

            for ( BwEdgeIterator e = _nodes[u]._in_edges_begin ; e != _nodes[u]._in_edges_end ; ++e )
            {
                assert( e != INVALID_BW_EDGE_ITERATOR );
                assert( e < _backward_edges.size() );
                assert( !_backward_edges[e].is_dummy() );
                assert( _backward_edges[e]._forward_edge != INVALID_FW_EDGE_ITERATOR );
                assert( !_forward_edges[_backward_edges[e]._forward_edge].is_dummy() );
                assert( _forward_edges[_backward_edges[e]._forward_edge]._backward_edge == e );
            }
        }

        return true;
    }

    bool has_parallel_edges() const
    {
        std::vector<BackwardEdge> edges;
        for ( NodeIterator u = nodes_begin() ; u != nodes_end() ; ++u )
            for ( BwEdgeIterator e = _nodes[u]._in_edges_begin ; e != _nodes[u]._in_edges_end ; ++e )
                edges.push_back(_backward_edges[e]);

        return
            util::has_duplicates
            (
                edges.begin(),
                edges.end(),
                [&] (const BackwardEdge& lhs, const BackwardEdge& rhs) -> bool // "less" operator
                {
                    assert( lhs._source < _nodes.size() );
                    assert( rhs._source < _nodes.size() );
                    assert( lhs._forward_edge < _forward_edges.size() );
                    assert( rhs._forward_edge < _forward_edges.size() );
                    assert( ! _forward_edges[lhs._forward_edge].is_dummy() );
                    assert( ! _forward_edges[rhs._forward_edge].is_dummy() );
                    assert( _forward_edges[lhs._forward_edge]._target < _nodes.size() );
                    assert( _forward_edges[rhs._forward_edge]._target < _nodes.size() );

                    if ( lhs._source < rhs._source )
                        return true;

                    if ( lhs._source == rhs._source )
                        if ( _forward_edges[lhs._forward_edge]._target < _forward_edges[rhs._forward_edge]._target )
                            return true;

                    return false;
                },
                [&] (const BackwardEdge& lhs, const BackwardEdge& rhs) -> bool // "equal_to" operator
                {
                    if ( lhs._source != rhs._source )
                        return false;

                    if ( _forward_edges[lhs._forward_edge]._target != _forward_edges[rhs._forward_edge]._target )
                        return false;

                    return true;
                }
       );
    }

    size_t get_n_nodes() const
    {
        return _nodes.size();
    }

    size_t get_n_edges() const
    {
        return _n_edges;
    }

    NodeIterator nodes_begin() const
    {
        return 0;
    }

    NodeIterator nodes_end() const
    {
        return _nodes.size();
    }

    FwEdgeIterator out_edges_begin(const NodeIterator src) const
    {
        assert( src < _nodes.size() );
        return _nodes[src]._out_edges_begin;
    }

    FwEdgeIterator out_edges_end(const NodeIterator src) const
    {
        assert( src < _nodes.size() );
        return _nodes[src]._out_edges_end;
    }

    NodeIterator get_target(const FwEdgeIterator e) const
    {
        assert( e < _forward_edges.size() );
        assert( ! _forward_edges[e].is_dummy() );
        return _forward_edges[e]._target;
    }

    BwEdgeIterator in_edges_begin(const NodeIterator tgt) const
    {
        assert( tgt < _nodes.size() );
        return _nodes[tgt]._in_edges_begin;
    }

    BwEdgeIterator in_edges_end(const NodeIterator tgt) const
    {
        assert( tgt < _nodes.size() );
        return _nodes[tgt]._in_edges_end;
    }

    NodeIterator get_source(const BwEdgeIterator e) const
    {
        assert( e < _backward_edges.size() );
        assert( ! _backward_edges[e].is_dummy() );
        return _backward_edges[e]._source;
    }

    const EdgeData& get_edge_data(FwEdgeIterator e) const
    {
        assert( e < _forward_edges.size() );
        assert( ! _forward_edges[e].is_dummy() );
        return _forward_edges[e]._data;
    }

    EdgeData& get_edge_data(FwEdgeIterator e)
    {
        return const_cast<EdgeData&>(static_cast<const DynamicGraph&>(*this).get_edge_data(e));
    }

    const EdgeData& get_edge_data(BwEdgeIterator e) const
    {
        assert( e < _backward_edges.size() );
        assert( ! _backward_edges[e].is_dummy() );
        FwEdgeIterator e_forward = _backward_edges[e]._forward_edge;

        assert( e_forward < _forward_edges.size() );
        assert( ! _forward_edges[e_forward].is_dummy() );
        return _forward_edges[e_forward]._data;
    }

    EdgeData& get_edge_data(BwEdgeIterator e)
    {
        return const_cast<EdgeData&>(static_cast<const DynamicGraph&>(*this).get_edge_data(e));
    }

    std::pair<FwEdgeIterator,BwEdgeIterator> get_edge_iterators(const NodeIterator src, const NodeIterator tgt) const
    {
        assert( src < _nodes.size() );
        assert( tgt < _nodes.size() );

        for ( BwEdgeIterator e = _nodes[tgt]._in_edges_begin ; e != _nodes[tgt]._in_edges_end ; ++e )
        {
            if ( _backward_edges[e]._source == src )
            {
                assert( _forward_edges[ _backward_edges[e]._forward_edge ]._target == tgt );
                return std::make_pair(_backward_edges[e]._forward_edge, e);
            }
        }

        return std::make_pair(INVALID_FW_EDGE_ITERATOR, INVALID_BW_EDGE_ITERATOR);
    }

    BwEdgeIterator get_backward_edge_iterator(const FwEdgeIterator& e_fw) const
    {
        assert( e_fw < _forward_edges.size() );
        assert( ! _forward_edges[e_fw].is_dummy() );
        assert( _forward_edges[e_fw]._backward_edge != INVALID_BW_EDGE_ITERATOR );
        assert( e_fw == _backward_edges[_forward_edges[e_fw]._backward_edge]._forward_edge );

        return _forward_edges[e_fw]._backward_edge;
    }

    FwEdgeIterator get_forward_edge_iterator(const BwEdgeIterator& e_bw) const
    {
        assert( e_bw < _backward_edges.size() );
        assert( ! _backward_edges[e_bw].is_dummy() );
        assert( _backward_edges[e_bw]._forward_edge != INVALID_FW_EDGE_ITERATOR );
        assert( e_bw == _forward_edges[_backward_edges[e_bw]._forward_edge]._backward_edge );

        return _backward_edges[e_bw]._forward_edge;
    }

    bool has_edge(const NodeIterator& src, const NodeIterator& tgt) const
    {
        for ( BwEdgeIterator e = _nodes[tgt]._in_edges_begin ; e != _nodes[tgt]._in_edges_end ; ++e )
            if ( _backward_edges[e]._source == src )
                return true;

        return false;
    }

    std::pair<FwEdgeIterator,BwEdgeIterator> insert_edge(const NodeIterator& src, const NodeIterator& tgt, EdgeData&& data)
    {
        assert( ! has_edge(src, tgt) );

        const FwEdgeIterator fw_insert_it =
                insert_edge_impl(_forward_edges, _nodes[src]._out_edges_begin, _nodes[src]._out_edges_end);

        const BwEdgeIterator bw_insert_it =
                insert_edge_impl(_backward_edges, _nodes[tgt]._in_edges_begin, _nodes[tgt]._in_edges_end);

        _forward_edges[fw_insert_it]._data = std::move(data);
        _forward_edges[fw_insert_it]._target = tgt;
        _forward_edges[fw_insert_it]._backward_edge = bw_insert_it;

        _backward_edges[bw_insert_it]._source = src;
        _backward_edges[bw_insert_it]._forward_edge = fw_insert_it;

        ++_n_edges;

        return std::make_pair(fw_insert_it, bw_insert_it);
    }

    void remove_all_edges_of_node(const NodeIterator& u)
    {

        assert( _n_edges >= _nodes[u]._in_edges_end - _nodes[u]._in_edges_begin );
        _n_edges -= _nodes[u]._in_edges_end - _nodes[u]._in_edges_begin;

        assert( _n_edges >= _nodes[u]._out_edges_end - _nodes[u]._out_edges_begin );
        _n_edges -= _nodes[u]._out_edges_end - _nodes[u]._out_edges_begin;

        //
        // remove outgoing edges and their corresponding backward edges
        //

        for ( FwEdgeIterator e_fw = _nodes[u]._out_edges_begin ; e_fw != _nodes[u]._out_edges_end ; ++e_fw )
        {
            assert( e_fw != INVALID_FW_EDGE_ITERATOR );
            assert( get_backward_edge_iterator(e_fw) != INVALID_BW_EDGE_ITERATOR );
            assert( _forward_edges[e_fw]._target != u );
            remove_in_edge(_forward_edges[e_fw]._target, get_backward_edge_iterator(e_fw));

            _forward_edges[e_fw].make_dummy();
        }

        _nodes[u]._out_edges_end = _nodes[u]._out_edges_begin;

        //
        // remove incoming edges and their corresponding forward edges
        //

        for ( BwEdgeIterator e_bw = _nodes[u]._in_edges_begin ; e_bw != _nodes[u]._in_edges_end ; ++e_bw )
        {
            assert( e_bw != INVALID_BW_EDGE_ITERATOR );
            assert( get_forward_edge_iterator(e_bw) != INVALID_FW_EDGE_ITERATOR );
            assert( _backward_edges[e_bw]._source != u );
            remove_out_edge(_backward_edges[e_bw]._source, get_forward_edge_iterator(e_bw));

            _backward_edges[e_bw].make_dummy();
        }

        _nodes[u]._in_edges_end = _nodes[u]._in_edges_begin;
    }

};

}

#endif /* KATCH_DYNAMIC_GRAPH_H_ */
