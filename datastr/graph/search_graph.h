/*
 * katch/datastr/graph/search_graph.h
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

#ifndef KATCH_SEARCH_GRAPH_H_
#define KATCH_SEARCH_GRAPH_H_

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include "datastr/base/ttf_wrapper.h"
#include "datastr/base/pwl_ttf.h"
#include "datastr/graph/basic.h"

namespace katch
{

class SearchGraph
{
private:

    static constexpr double PERIOD = 864000.0;
    static constexpr int BUCKET_SIZE = 8;

    using TTFImpl = PwlTTF<static_cast<int>(PERIOD), BUCKET_SIZE>;

public:

    using TTFRef = TTFReference<TTFImpl>;
    using TTF = TTFWrapper<TTFImpl>;

private:

    struct Node
    {
        EdgeIterator _first_edge;
        Node() : _first_edge(INVALID_EDGE_ITERATOR) {}
    };

    static constexpr NodeIterator UNINITIALIZED_NODE_IT = (1 << 28) - 1;

    struct Edge
    {
        NodeIterator _other_node : 28;
        bool _upward : 1;
        bool _downward : 1;
        bool _is_exact : 1;
        bool _is_constant : 1;
        uint32_t _first_middle_node;

        union
        {
            double _constant_value;
            const TTFImpl* _ttf_impl_ptr;
            struct { uint32_t _lower; uint32_t _upper; };
        };

        Edge(const Edge&) = delete;
        Edge& operator= (const Edge&) = delete;

        Edge(Edge&& edge)
        {
            _other_node = edge._other_node;
            _upward = edge._upward;
            _downward = edge._downward;
            _is_exact = edge._is_exact;
            _is_constant = edge._is_constant;
            _first_middle_node = edge._first_middle_node;

            if ( _is_constant && _is_exact )
            {
                _constant_value = edge._constant_value;
            }
            else if ( _is_constant && ! _is_exact )
            {
                _lower = edge._lower;
                _upper = edge._upper;
            }
            else
            {
                assert( ! _is_constant );

                _ttf_impl_ptr = edge._ttf_impl_ptr;
                edge._ttf_impl_ptr = nullptr;
            }
        }

        Edge& operator= (Edge&& edge)
        {
            if ( this != &edge )
            {
                if ( ! _is_constant )
                    if ( _ttf_impl_ptr )
                        delete _ttf_impl_ptr;

                _other_node = edge._other_node;
                _upward = edge._upward;
                _downward = edge._downward;
                _is_exact = edge._is_exact;
                _is_constant = edge._is_constant;
                _first_middle_node = edge._first_middle_node;

                if ( _is_constant && _is_exact )
                {
                    _constant_value = edge._constant_value;
                }
                else if ( _is_constant && ! _is_exact )
                {
                    _lower = edge._lower;
                    _upper = edge._upper;
                }
                else
                {
                    assert( ! _is_constant );

                    _ttf_impl_ptr = edge._ttf_impl_ptr;
                    edge._ttf_impl_ptr = nullptr;
                }
            }

            return *this;
        }

        Edge()
        : _other_node(UNINITIALIZED_NODE_IT),
          _upward(false),
          _downward(false),
          _is_exact(true),
          _is_constant(false),
          _first_middle_node(0),
          _ttf_impl_ptr(nullptr)
        {}

        ~Edge()
        {
            if ( ! _is_constant )
                if ( _ttf_impl_ptr )
                    delete _ttf_impl_ptr;

            _is_constant = false;
            _ttf_impl_ptr = nullptr;
        }
    };

    std::vector<Node> _nodes;
    std::vector<Edge> _edges;
    std::vector<MiddleNodeDescr> _middle_node_data;

    double _lower_bound_factor;

    template <typename EdgeInfo>
    static void bidirectionalize(std::vector<EdgeInfo>& edge_list)
    {
        std::sort(
            edge_list.begin(), edge_list.end(),
                [](const EdgeInfo& lhs, const EdgeInfo& rhs) -> bool
                {
                    if ( lhs.get_source() < rhs.get_source() )
                        return true;

                    if ( lhs.get_source() == rhs.get_source() )
                    {
                        if ( lhs.get_target() < rhs.get_target() )
                            return true;

                        if ( lhs.get_target() == rhs.get_target() )
                            return lhs.get_ttf().is_constant() && ! rhs.get_ttf().is_constant();
                    }

                    return false;
                }
        );

        for ( size_t i = 0 ; i + 1 < edge_list.size() ; ++i )
        {
            if
            (
                    edge_list[i].get_source() == edge_list[i+1].get_source() &&
                    edge_list[i].get_target() == edge_list[i+1].get_target() &&
                    edge_list[i].get_forward() != edge_list[i+1].get_forward() &&
                    edge_list[i].get_backward() != edge_list[i+1].get_backward() &&
                    edge_list[i].get_ttf().is_constant() && edge_list[i+1].get_ttf().is_constant() &&
                    edge_list[i].get_ttf().get_constant_value() == edge_list[i+1].get_ttf().get_constant_value() &&
                    edge_list[i].get_middle_node_data() == edge_list[i+1].get_middle_node_data()
            )
            {
                edge_list[i].set_forward(true);
                edge_list[i].set_backward(true);
                edge_list[i+1].set_forward(true);
                edge_list[i+1].set_backward(true);
            }
        }

        auto new_end =
            std::unique(
                edge_list.begin(), edge_list.end(),
                    [](const EdgeInfo& lhs, const EdgeInfo& rhs) -> bool
                    {
                        if ( lhs.get_source() != rhs.get_source() )
                            return false;

                        if ( lhs.get_target() != rhs.get_target() )
                            return false;

                        if ( ! (lhs.get_forward() && rhs.get_forward() && lhs.get_backward() && rhs.get_backward() ) )
                            return false;

                        return true;
                    }
            );

        edge_list.resize( std::distance(edge_list.begin(), new_end) );
    }

public:

    template <typename EdgeInfo>
    SearchGraph(std::vector<EdgeInfo>&& edge_list, const double epsilon = 0.0)
    : _nodes(), _edges(), _middle_node_data(),
      _lower_bound_factor( 1.0 / (1.0 + epsilon) )
    {
        bidirectionalize(edge_list);

        size_t max_node_it = 0;
        std::for_each( edge_list.begin(), edge_list.end(),
                        [&max_node_it](const EdgeInfo& edge_info) -> void
                        {
                            max_node_it = std::max(max_node_it, edge_info.get_source());
                            max_node_it = std::max(max_node_it, edge_info.get_target());
                        } );

        _nodes.assign(max_node_it + 2, Node());

        std::sort( edge_list.begin(), edge_list.end(),
                        [](const EdgeInfo& lhs, const EdgeInfo& rhs) -> bool
                                { return lhs.get_source() < rhs.get_source(); } );

        assert( _nodes[0]._first_edge == INVALID_EDGE_ITERATOR );
        _nodes[0]._first_edge = 0;

        NodeIterator current_source = 0;

        for ( auto it = edge_list.begin() ; it != edge_list.end() ; ++it )
        {
            assert( it->get_source() + 1 < _nodes.size() );

            while ( current_source <= it->get_source() )
                _nodes[++current_source]._first_edge = _edges.size();

            _edges.push_back(std::move(Edge()));
            _edges.back()._other_node = it->get_target();
            _edges.back()._upward = it->get_forward();
            _edges.back()._downward = it->get_backward();
            _edges.back()._is_exact = true;

            _edges.back()._first_middle_node = _middle_node_data.size();

            _middle_node_data.insert(
                    _middle_node_data.end(), it->get_middle_node_data().begin(), it->get_middle_node_data().end());

            if ( it->get_ttf().is_constant() )
            {
                _edges.back()._is_constant = true;
                _edges.back()._constant_value = it->get_ttf().get_constant_value();
            }
            else
            {
                assert( it->get_ttf().get_ttf_impl_ptr() );

                _edges.back()._is_constant = false;
                _edges.back()._ttf_impl_ptr = it->get_ttf().get_ttf_impl_ptr();
                it->get_ttf().release();
            }

            ++_nodes[current_source]._first_edge;
        }

        while ( current_source <= max_node_it )
            _nodes[++current_source]._first_edge = _edges.size();

        assert( _nodes.back()._first_edge == _edges.size() );

        _edges.emplace_back();
        _edges.back()._first_middle_node = _middle_node_data.size();

        edge_list.clear();
    }

    size_t get_n_nodes() const noexcept { return _nodes.size() - 1; }
    size_t get_n_edges() const noexcept { return _edges.size() - 1; }

    NodeIterator nodes_begin() const noexcept { return 0; }
    NodeIterator nodes_end() const noexcept { return _nodes.size() - 1; }

    EdgeIterator edges_begin(const NodeIterator u) const noexcept
    {
        assert( u + 1 < _nodes.size() );
        return _nodes[u]._first_edge;
    }

    EdgeIterator edges_end(const NodeIterator u) const noexcept
    {
        assert( u + 1 < _nodes.size() );
        return _nodes[u+1]._first_edge;
    }

    NodeIterator get_other_node(const EdgeIterator e) const noexcept
    {
        assert( e + 1 < _edges.size() );
        return _edges[e]._other_node;
    }

    bool is_directed_upward(const EdgeIterator e) const noexcept
    {
        assert( e < _edges.size() );
        return _edges[e]._upward;
    }

    bool is_directed_downward(const EdgeIterator e) const noexcept
    {
        assert( e + 1 < _edges.size() );
        return _edges[e]._downward;
    }

    double get_upper(const EdgeIterator e) const noexcept
    {
        assert( e + 1 < _edges.size() );
        if ( _edges[e]._is_exact )
        {
            if ( _edges[e]._is_constant )
                return _edges[e]._constant_value;

            assert( _edges[e]._ttf_impl_ptr );
            return _edges[e]._ttf_impl_ptr->get_max();
        }

        if ( _edges[e]._is_constant )
            return double(_edges[e]._upper);

        return _edges[e]._ttf_impl_ptr->get_max();
    }

    double get_lower(const EdgeIterator e) const noexcept
    {
        assert( e + 1 < _edges.size() );
        if ( _edges[e]._is_exact )
        {
            if ( _edges[e]._is_constant )
                return _edges[e]._constant_value;

            assert( _edges[e]._ttf_impl_ptr );
            return _edges[e]._ttf_impl_ptr->get_min();
        }

        if ( _edges[e]._is_constant )
            return double(_edges[e]._lower);

        return _lower_bound_factor * _edges[e]._ttf_impl_ptr->get_min();
    }

    TTFRef get_ttf(const EdgeIterator e) const noexcept
    {
        assert( e + 1 < _edges.size() );
        assert( _edges[e]._is_exact );

        if ( _edges[e]._is_constant )
            return TTFRef(_edges[e]._constant_value);

        assert( _edges[e]._ttf_impl_ptr );
        return TTFRef(_edges[e]._ttf_impl_ptr);
    }

    TTFRef get_upper_bound(const EdgeIterator e) const noexcept
    {
        assert( e + 1 < _edges.size() );
        assert( ! _edges[e]._is_exact );

        if ( _edges[e]._is_constant )
            return TTFRef(double(_edges[e]._upper));

        assert( _edges[e]._ttf_impl_ptr );
        return TTFRef(_edges[e]._ttf_impl_ptr);
    }

    TTF get_lower_bound(const EdgeIterator e) const noexcept
    {
        assert( e + 1 < _edges.size() );
        assert( ! _edges[e]._is_exact );

        if ( _edges[e]._is_constant )
            return TTF(double(_edges[e]._lower));

        assert( _edges[e]._ttf_impl_ptr );
        return TTF(_edges[e]._ttf_impl_ptr->scale(_lower_bound_factor));
    }

    size_t get_n_middle_nodes(const EdgeIterator e) const noexcept
    {
        assert( e + 1 < _edges.size() );
        return _edges[e+1]._first_middle_node - _edges[e]._first_middle_node;
    }

    MiddleNodeDescr get_middle_node_descr(const EdgeIterator e, const size_t index) const noexcept
    {
        assert( e + 1 < _edges.size() );

        if ( get_n_middle_nodes(e) == 0 )
            return MiddleNodeDescr(0.0, INVALID_NODE_ITERATOR);

        assert( index < get_n_middle_nodes(e) );
        const size_t position = _edges[e]._first_middle_node + index;

        assert( position < _edges[e+1]._first_middle_node );
        return _middle_node_data[position];
    }

    NodeIterator get_middle_node(const EdgeIterator e, const double time) const noexcept
    {
        assert( e + 1 < _edges.size() );

        if ( get_n_middle_nodes(e) == 0 )
            return INVALID_NODE_ITERATOR;

        if ( get_n_middle_nodes(e) == 1 )
            return _middle_node_data[_edges[e]._first_middle_node]._middle_node;

        const double x = util::modulo(time, PERIOD);

        assert( _edges[e]._first_middle_node <= _middle_node_data.size() );
        assert( _edges[e+1]._first_middle_node <= _middle_node_data.size() );
        assert( _edges[e]._first_middle_node <= _edges[e+1]._first_middle_node );

        uint32_t index;
        for ( index = _edges[e]._first_middle_node ; index != _edges[e+1]._first_middle_node ; ++index )
            if ( x < _middle_node_data[index]._time )
                break;

        if ( index == _edges[e]._first_middle_node )
        {
            assert( _edges[e+1]._first_middle_node > 0 );
            return _middle_node_data[_edges[e+1]._first_middle_node - 1]._middle_node;
        }

        if ( index + 1 == _edges[e+1]._first_middle_node )
        {
            assert( _edges[e+1]._first_middle_node > _edges[e]._first_middle_node );

            if ( le(_middle_node_data[_edges[e+1]._first_middle_node - 1]._time, x) )
                return _middle_node_data[_edges[e+1]._first_middle_node - 1]._middle_node;

            assert( x < _middle_node_data[_edges[e+1]._first_middle_node - 1]._time );
        }

        assert( index > 0 );
        assert( index - 1 >= _edges[e]._first_middle_node );
        assert( index - 1 < _edges[e+1]._first_middle_node );

        assert( le(_middle_node_data[index - 1]._time, x) );
        assert( KATCH_IMPLIES(index < _edges[e+1]._first_middle_node, le( x, _middle_node_data[index]._time)) );

        return _middle_node_data[index - 1]._middle_node;
    }
};

}



#endif /* KATCH_SEARCH_GRAPH_H_ */
