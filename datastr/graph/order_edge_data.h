/*
 * katch/datastr/graph/order_edge_data.h
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

#ifndef KATCH_ORDER_EDGE_DATA_H_
#define KATCH_ORDER_EDGE_DATA_H_

#include <assert.h>
#include <limits>
#include <memory>
#include <utility>

#include "datastr/base/pwl_ttf.h"
#include "datastr/base/ttf_wrapper.h"
#include "datastr/graph/basic.h"
#include "util/misc.h"

namespace katch
{

class OrderEdgeData
{

private:

    static constexpr double PERIOD = 864000.0;
    static constexpr int BUCKET_SIZE = 8;

public:

    using TTFImpl = PwlTTF<static_cast<int>(PERIOD), BUCKET_SIZE>;
    using TTFRef = TTFReference<TTFImpl>;
    using TTF = TTFWrapper<TTFImpl>;
    using UndercutDescriptor = TTF::UndercutDescriptor;

private:

    // number of original edges represented by this (shortcut) edge
    unsigned int _n_original_edges;

    // size of of shortcut descriptor
    unsigned int _shortcut_descr_size: 31;

    // tells whether this edge is directed forward
    bool _is_constant : 1;

    union
    {
        // this edges TTF if it is constant
        double _constant_value;

        // this edges TTF if it is not constant
        const TTFImpl* _ttf_impl_ptr;
    };

    union
    {
        // the middle node if it is unique
        NodeIterator _unique_middle_node;

        // shortcut data
        MiddleNodeDescr* _middle_node_data;
    };


public:

    OrderEdgeData(const OrderEdgeData&) = delete;
    OrderEdgeData& operator= (const OrderEdgeData&) = delete;

    OrderEdgeData(OrderEdgeData&& order_edge_data)
    : _n_original_edges(order_edge_data._n_original_edges),
      _shortcut_descr_size(order_edge_data._shortcut_descr_size),
      _is_constant(order_edge_data._is_constant)
    {
        if ( _is_constant )
            _constant_value = order_edge_data._constant_value;
        else
            _ttf_impl_ptr = order_edge_data._ttf_impl_ptr;

        if ( _shortcut_descr_size == 1 )
            _unique_middle_node = order_edge_data._unique_middle_node;
        else
            _middle_node_data = order_edge_data._middle_node_data;

        order_edge_data._is_constant = false;
        order_edge_data._ttf_impl_ptr = nullptr;

        order_edge_data._shortcut_descr_size = 0;
        order_edge_data._middle_node_data = nullptr;
    }

    OrderEdgeData& operator= (OrderEdgeData&& order_edge_data)
    {
        if ( this != &order_edge_data )
        {
            if ( ! _is_constant )
                if ( _ttf_impl_ptr )
                    delete _ttf_impl_ptr;

            if ( _shortcut_descr_size > 1 )
                if ( _middle_node_data )
                    free( _middle_node_data );

            _n_original_edges = order_edge_data._n_original_edges;
            _shortcut_descr_size = order_edge_data._shortcut_descr_size;
            _is_constant = order_edge_data._is_constant;

            if ( _is_constant )
                _constant_value = order_edge_data._constant_value;
            else
                _ttf_impl_ptr = order_edge_data._ttf_impl_ptr;

            if ( _shortcut_descr_size == 1 )
                _unique_middle_node = order_edge_data._unique_middle_node;
            else
                _middle_node_data = order_edge_data._middle_node_data;

            order_edge_data._is_constant = false;
            order_edge_data._ttf_impl_ptr = nullptr;

            order_edge_data._shortcut_descr_size = 0;
            order_edge_data._middle_node_data = nullptr;
        }

        return *this;
    }

    OrderEdgeData()
    : _n_original_edges(1),
      _shortcut_descr_size(0),
      _is_constant(false),
      _ttf_impl_ptr(nullptr),
      _middle_node_data(nullptr)
    {}

    ~OrderEdgeData()
    {
        if ( ! _is_constant )
            if ( _ttf_impl_ptr )
                delete _ttf_impl_ptr;

        _is_constant = false;
        _ttf_impl_ptr = nullptr;

        if ( _shortcut_descr_size > 1 )
            if ( _middle_node_data )
                free( _middle_node_data );

        _shortcut_descr_size = 0;
        _middle_node_data = nullptr;
    }

    template <typename TTFType>
    void set_ttf(TTFType&& ttf)
    {
        if ( ! _is_constant )
            if ( _ttf_impl_ptr )
                delete _ttf_impl_ptr;

        if ( ttf.is_constant() )
        {
            _is_constant = true;
            _constant_value = ttf.get_constant_value();
        }
        else
        {
            _is_constant = false;
            _ttf_impl_ptr = ttf.get_ttf_impl_ptr();
        }

        ttf.release();
    }

    template <typename EdgeInfoT>
    void set_up(EdgeInfoT&& edge_info) noexcept
    {
        set_ttf( std::move(edge_info.get_ttf()) );
    }

    TTFRef get_ttf() const noexcept
    {
        if ( _is_constant )
            return TTFRef(_constant_value);

        return TTFRef(_ttf_impl_ptr);
    }

    void set_n_original_edges(const unsigned int n_original_edges) noexcept
    {
        _n_original_edges = n_original_edges;
    }

    unsigned int get_n_original_edges() const noexcept
    {
        return _n_original_edges;
    }

    bool never_represents_a_shortcut() const noexcept
    {
        if ( _shortcut_descr_size == 0 )
            return true;

        return false;
    }

    size_t get_n_middle_nodes() const
    {
        return _shortcut_descr_size;
    }

    MiddleNodeDescr get_middle_node_descr(int index) const
    {
        const int orig_index = index;

        if ( _shortcut_descr_size == 0 )
        {
            double time = 0.0;
            if ( index == 1 ) time = PERIOD;
            else if ( index == -1 ) time = -PERIOD;
            else assert( index == 0 );

            return MiddleNodeDescr(time, INVALID_NODE_ITERATOR);
        }

        if ( _shortcut_descr_size == 1 )
        {
            double time = 0.0;
            if ( index == 1 ) time = PERIOD;
            else if ( index == -1 ) time = -PERIOD;
            else assert( index == 0 );

            return MiddleNodeDescr(time, _unique_middle_node);
        }

        double offset = 0.0;

        if ( index >= -((int)_shortcut_descr_size) && index < (int) 0 )
        {
            index += (int) _shortcut_descr_size;
            offset =  -PERIOD;
        }
        else if ( index >= (int) _shortcut_descr_size && index < (int) 2*_shortcut_descr_size )
        {
            index -= (int) _shortcut_descr_size;
            offset =  PERIOD;
        }

        return MiddleNodeDescr(_middle_node_data[index]._time + offset, _middle_node_data[index]._middle_node);
    }

    NodeIterator get_middle_node(const double time) const noexcept
    {
        if ( _shortcut_descr_size == 0 )
            return INVALID_NODE_ITERATOR;

        if ( _shortcut_descr_size == 1 )
            return _unique_middle_node;

        const double x = util::modulo(time, PERIOD);

        uint32_t index;
        for ( index = 0 ; index != _shortcut_descr_size ; ++index )
            if ( x < _middle_node_data[index]._time )
                break;

        if ( index == 0 )
        {
            assert( _shortcut_descr_size > 0 );
            return _middle_node_data[_shortcut_descr_size - 1]._middle_node;
        }

        if ( index + 1 == _shortcut_descr_size )
        {
            if ( le(_middle_node_data[_shortcut_descr_size - 1]._time, x) )
                return _middle_node_data[_shortcut_descr_size - 1]._middle_node;

            assert( x < _middle_node_data[_shortcut_descr_size - 1]._time );
        }

        assert( index > 0 );
        assert( index - 1 < _shortcut_descr_size );

        assert( le(_middle_node_data[index - 1]._time, x) );
        assert( KATCH_IMPLIES(index <  _shortcut_descr_size, le( x, _middle_node_data[index]._time)) );

        return _middle_node_data[index - 1]._middle_node;
    }

    void set_original_edge()
    {
        if ( _shortcut_descr_size > 1 )
            if ( _middle_node_data )
                free( _middle_node_data );

        _shortcut_descr_size = 0;
        _middle_node_data = nullptr;
    }

    void set_unique_middle_node(const NodeIterator& middle_node)
    {
        if ( _shortcut_descr_size > 1 )
            if ( _middle_node_data )
                free( _middle_node_data );

        _shortcut_descr_size = 1;
        _unique_middle_node = middle_node;
    }

    template <typename Iterator>
    void set_middle_node_data(Iterator begin, Iterator end)
    {
        if ( end == begin )
        {
            set_original_edge();
            return;
        }

        if ( end - begin == 1 )
        {
            set_unique_middle_node(begin->_middle_node);
            return;
        }

        if (
                std::all_of(begin, end,
                        [](const MiddleNodeDescr& descr){ return descr._middle_node == INVALID_NODE_ITERATOR; })
        )
        {
            set_original_edge();
            return;
        }

        assert( end - begin >= 2 );

        if ( _shortcut_descr_size > 1 )
            if ( _middle_node_data )
                free( _middle_node_data );

        _shortcut_descr_size = end - begin;
        _middle_node_data = (MiddleNodeDescr*) malloc( _shortcut_descr_size * sizeof(*_middle_node_data) );

        for ( auto it = begin ; it != end ; ++it )
        {
            assert( it - begin < _shortcut_descr_size );
            _middle_node_data[it - begin] = *it;
        }
    }

};

}

#endif /* KATCH_ORDER_EDGE_DATA_H_ */
