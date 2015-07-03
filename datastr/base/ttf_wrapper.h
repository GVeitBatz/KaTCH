/*
 * katch/datastr/base/ttf_wrapper.h
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

#ifndef KATCH_TTF_WRAPPER_H_
#define KATCH_TTF_WRAPPER_H_

#include <assert.h>
#include <limits>
#include <utility>
#include <deque>

#include "datastr/base/double.h"
#include "util/misc.h"

namespace katch
{

template <typename ResultTTF, typename LhsTTF, typename RhsTTF>
inline ResultTTF link_impl(const LhsTTF& f, const RhsTTF& g) noexcept
{
    if ( f.is_constant() && g.is_constant() )
        return ResultTTF( f.get_constant_value() + g.get_constant_value() );

    if ( f.is_constant() && ! g.is_constant() )
        return ResultTTF( link(f.get_constant_value(), g.get_ttf_impl_ptr()) );

    if ( ! f.is_constant() && g.is_constant() )
        return ResultTTF( link(f.get_ttf_impl_ptr(), g.get_constant_value()) );

    assert( ! f.is_constant() && ! g.is_constant() );
    return ResultTTF( link(f.get_ttf_impl_ptr(), g.get_ttf_impl_ptr()) );
}


template <typename ResultTTF, typename LhsTTF, typename RhsTTF, typename UndercutDescr>
inline ResultTTF merge_impl(const LhsTTF& f, const RhsTTF& g, UndercutDescr& descr) noexcept
{
    if ( f.is_constant() && g.is_constant() )
    {
        descr._f_undercuts_strictly = lt(f.get_constant_value(), g.get_constant_value());
        descr._g_undercuts_strictly = lt(g.get_constant_value(), f.get_constant_value());

        return ResultTTF( std::min(f.get_constant_value(), g.get_constant_value()) );
    }

    if ( f.is_constant() && ! g.is_constant() )
    {
        if ( lt(f.get_constant_value(), g.get_min()) )
        {
            descr._f_undercuts_strictly = true;
            descr._g_undercuts_strictly = false;
            return ResultTTF(f.get_constant_value());
        }
        else if ( lt(g.get_max(), f.get_constant_value()) )
        {
            descr._f_undercuts_strictly = false;
            descr._g_undercuts_strictly = true;
            return g.clone();
        }

        return ResultTTF( merge(f.get_constant_value(), g.get_ttf_impl_ptr(), descr) );
    }

    if ( ! f.is_constant() && g.is_constant() )
    {
        if ( lt(g.get_constant_value(), f.get_min()) )
        {
            descr._g_undercuts_strictly = true;
            descr._f_undercuts_strictly = false;
            return ResultTTF(g.get_constant_value());
        }
        else if ( lt(f.get_max(), g.get_constant_value()) )
        {
            descr._f_undercuts_strictly = true;
            descr._g_undercuts_strictly = false;
            return f.clone();
        }

        return ResultTTF( merge(f.get_ttf_impl_ptr(), g.get_constant_value(), descr) );
    }

    assert( ! f.is_constant() && ! g.is_constant() );

    if ( lt(f.get_max(), g.get_min()) )
    {
        descr._f_undercuts_strictly = true;
        descr._g_undercuts_strictly = false;
        return f.clone();
    }
    else if ( lt(g.get_max(), f.get_min()) )
    {
        descr._f_undercuts_strictly = false;
        descr._g_undercuts_strictly = true;
        return g.clone();
    }

    return ResultTTF( merge(f.get_ttf_impl_ptr(), g.get_ttf_impl_ptr(), descr) );
}


template <typename TTFImplT>
class TTFWrapper;


template <typename TTFImplT>
class TTFReference
{

    friend class TTFWrapper<TTFImplT>;

public:

    using TTFImpl = TTFImplT;
    using UndercutDescriptor = typename TTFImpl::UndercutDescriptor;

    static constexpr double PERIOD() { return TTFImpl::PERIOD(); }

private:

    bool _is_constant;

    union
    {
        double _constant_value;
        const TTFImpl* _ttf_impl_ptr;
    };


public:

    TTFReference(TTFReference&& wrapper)
    : _is_constant(wrapper._is_constant)
    {
        if ( wrapper._is_constant )
           _constant_value = wrapper._constant_value;
        else
           _ttf_impl_ptr = wrapper._ttf_impl_ptr;

        #ifndef DNDEBUG
            wrapper._is_constant = false;
            wrapper._ttf_impl_ptr = nullptr;
        #endif
    }

    TTFReference& operator= (TTFReference&& other)
    {
        if ( this != &other )
        {
            if ( other._is_constant )
                _constant_value = other._constant_value;
            else
                _ttf_impl_ptr = other._ttf_impl_ptr;

            _is_constant = other._is_constant;

            #ifndef DNDEBUG
                other._is_constant = false;
                other._ttf_impl_ptr = nullptr;
            #endif
        }

        return *this;
    }

    TTFReference(const TTFReference& ttf_ref)
    : _is_constant(ttf_ref._is_constant)
    {
        if ( ttf_ref._is_constant )
            _constant_value = ttf_ref._constant_value;
        else
            _ttf_impl_ptr = ttf_ref._ttf_impl_ptr;
    }

    TTFReference& operator= (const TTFReference& other)
    {
        if ( this != &other )
        {
            if ( other._is_constant )
                _constant_value = other._constant_value;
            else
                _ttf_impl_ptr = other._ttf_impl_ptr;

            _is_constant = other._is_constant;
        }

        return *this;
    }

    TTFReference(const double& c)
    : _is_constant(true), _constant_value(c)
    {}

    TTFReference(const TTFImpl* ttf_impl_ptr)
    : _is_constant(false), _ttf_impl_ptr(ttf_impl_ptr)
    {}

    TTFReference()
    : TTFReference( std::numeric_limits<double>::max() )
    {}

    const Point operator[] (const int& index) const noexcept
    {
        if ( _is_constant )
        {
            if ( index == 0 )
                return Point(0.0, _constant_value);

            if ( index == 1 )
                return Point(PERIOD(), _constant_value);

            if ( index == -1 )
                return Point(-PERIOD(), _constant_value);

            return Point( PERIOD() * double(index), _constant_value );
        }

        return _ttf_impl_ptr->operator[] (index);
    }

    bool is_constant() const noexcept
    {
        return _is_constant;
    }

    double eval(const double& x) const noexcept
    {
        if ( _is_constant ) return _constant_value;

        assert( _ttf_impl_ptr );
        return _ttf_impl_ptr->eval(x);
    }

    size_t size() const noexcept
    {
        if ( _is_constant ) return 0;

        assert( _ttf_impl_ptr );
        return _ttf_impl_ptr->size();
    }

    double get_constant_value() const
    {
        assert( _is_constant );
        return _constant_value;
    }

    const TTFImpl* get_ttf_impl_ptr() const
    {
        assert( ! _is_constant );
        return _ttf_impl_ptr;
    }

    double get_min() const noexcept
    {
        if ( _is_constant ) return _constant_value;

        assert( _ttf_impl_ptr );
        return _ttf_impl_ptr->get_min();
    }

    double get_max() const noexcept
    {
        if ( _is_constant ) return _constant_value;

        assert( _ttf_impl_ptr );
        return _ttf_impl_ptr->get_max();
    }

    std::pair<double,double> get_min_max() const noexcept
    {
        if ( _is_constant ) return _constant_value;

        assert( _ttf_impl_ptr );
        return std::make_pair(_ttf_impl_ptr->get_min(), _ttf_impl_ptr->get_max());
    }

    TTFWrapper<TTFImpl> add(const double& c) const noexcept
    {
        if ( _is_constant )
            return TTFWrapper<TTFImpl>(_constant_value + c);

        return TTFWrapper<TTFImpl>(_ttf_impl_ptr->add(c));
    }

    TTFWrapper<TTFImpl> clone() const noexcept
    {
        if ( _is_constant )
            return TTFWrapper<TTFImpl>( _constant_value );

        return TTFWrapper<TTFImpl>(new TTFImpl(*_ttf_impl_ptr));
    }

    template <typename ResultTTF, typename LhsTTF, typename RhsTTF>
    friend ResultTTF link_impl(const LhsTTF&, const RhsTTF&) noexcept;

    template <typename ResultTTF, typename LhsTTF, typename RhsTTF, typename UndercutDescr>
    friend ResultTTF merge_impl(const LhsTTF&, const RhsTTF&, UndercutDescr&) noexcept;

    friend std::ostream& operator<<(std::ostream& os, const TTFReference& f)
    {
        if ( f._is_constant )
            os << "{ const = " << f._constant_value << " }";
        else
            os << *f._ttf_impl_ptr;

        return os;
    }

    friend bool operator== (const TTFReference& lhs, const TTFReference& rhs)
    {
        if ( lhs._is_constant != rhs._is_constant ) return false;

        if ( lhs._is_constant )
            return lhs._constant_value == rhs._constant_value;

        assert( lhs._ttf_impl_ptr != rhs._ttf_impl_ptr );
        return false;
    }

    friend bool operator!= (const TTFReference& lhs, const TTFReference& rhs)
    {
        return ! (lhs == rhs);
    }
};





template <typename TTFImplT>
class TTFWrapper
{
public:

    using TTFImpl = TTFImplT;
    using UndercutDescriptor = typename TTFImpl::UndercutDescriptor;

    static constexpr double PERIOD() { return TTFImpl::PERIOD(); }
    static const TTFWrapper INFTY;

private:

    bool dbg_is_a_hidden_constant_function() const
    {
        if ( _ttf._is_constant )
            return false;

        if ( _ttf._ttf_impl_ptr->get_n_points() == 1 )
            return true;

        if ( _ttf._ttf_impl_ptr->get_n_points() == 2 )
            if ( eq( (*_ttf._ttf_impl_ptr)[0].y, (*_ttf._ttf_impl_ptr)[1].y) )
                return true;

        return false;
    }

    TTFReference<TTFImpl> _ttf;

public:

    TTFWrapper(const TTFWrapper&) = delete;
    TTFWrapper& operator= (const TTFWrapper&) = delete;

    TTFWrapper(const double& c) : _ttf(c) {}
    TTFWrapper(const TTFImpl* ttf_ptr) : _ttf(ttf_ptr) { assert( ttf_ptr != nullptr ); }

    TTFWrapper() noexcept
    : TTFWrapper( std::numeric_limits<double>::max() )
    {
        assert( *this == INFTY );
    }

    template <typename PointIterator>
    TTFWrapper(PointIterator points_begin, PointIterator points_end) noexcept
    : _ttf(new TTFImpl(points_begin, points_end))
    {
        assert( points_end >= points_begin );

        if ( _ttf.size() == 1 )
        {
            double constant_value = _ttf[0].y;

            delete _ttf._ttf_impl_ptr;

            _ttf._is_constant = true;
            _ttf._constant_value = constant_value;
        }

        assert( ! dbg_is_a_hidden_constant_function() );
    }

    TTFWrapper(TTFWrapper&& other) noexcept
    {
        if ( other._ttf._is_constant )
            _ttf._constant_value = other._ttf._constant_value;
        else
            _ttf._ttf_impl_ptr = other._ttf._ttf_impl_ptr;

        _ttf._is_constant = other._ttf._is_constant;

        other._ttf._is_constant = false;
        other._ttf._ttf_impl_ptr = nullptr;
    }

    TTFWrapper& operator= (TTFWrapper&& other) noexcept
    {
        if ( this != &other )
        {
            if ( ! _ttf._is_constant && _ttf._ttf_impl_ptr )
                delete _ttf._ttf_impl_ptr;

            if ( other._ttf._is_constant )
                _ttf._constant_value = other._ttf._constant_value;
            else
                _ttf._ttf_impl_ptr = other._ttf._ttf_impl_ptr;

            _ttf._is_constant = other._ttf._is_constant;

            other._ttf._is_constant = false;
            other._ttf._ttf_impl_ptr = nullptr;
        }

        return *this;
    }

    ~TTFWrapper()
    {
        if ( ! _ttf._is_constant )
            if ( _ttf._ttf_impl_ptr )
                delete _ttf._ttf_impl_ptr;

        _ttf._ttf_impl_ptr = nullptr;

    }

    void release() { _ttf._ttf_impl_ptr = nullptr; }

    double get_constant_value() const { return _ttf.get_constant_value(); }
    const TTFImpl* get_ttf_impl_ptr() const { return _ttf.get_ttf_impl_ptr(); }

    const Point operator[] (const int& index) const noexcept { return _ttf[index]; }
    bool is_constant() const noexcept { return _ttf.is_constant(); }
    double eval(const double& x) const noexcept { return _ttf.eval(x); }
    size_t size() const noexcept { return _ttf.size(); }
    double get_min() const noexcept { return _ttf.get_min(); }
    double get_max() const noexcept { return _ttf.get_max(); }
    std::pair<double,double> get_min_max() const noexcept { return _ttf.get_min_max(); }
    TTFWrapper add(const double& c) const noexcept { return _ttf.add(c); }
    TTFWrapper clone() const noexcept { return _ttf.clone(); }

    friend std::ostream& operator<<(std::ostream& os, const TTFWrapper& f) { return os << f._ttf; }
    friend bool operator== (const TTFWrapper& lhs, const TTFWrapper& rhs) { return lhs._ttf == rhs._ttf; }
    friend bool operator!= (const TTFWrapper& lhs, const TTFWrapper& rhs) { return ! (lhs == rhs); }

    template <typename ResultTTF, typename LhsTTF, typename RhsTTF>
    friend ResultTTF link_impl(const LhsTTF&, const RhsTTF&) noexcept;

    template <typename ResultTTF, typename LhsTTF, typename RhsTTF, typename UndercutDescr>
    friend ResultTTF merge_impl(const LhsTTF&, const RhsTTF&, UndercutDescr&) noexcept;

};

template <typename TTFImplT>
const TTFWrapper<TTFImplT> TTFWrapper<TTFImplT>::INFTY(std::numeric_limits<double>::max());

template <typename LhsTTF, typename RhsTTF>
TTFWrapper<typename LhsTTF::TTFImpl> link(const LhsTTF& f, const RhsTTF& g)
{
    return link_impl<TTFWrapper<typename LhsTTF::TTFImpl>, LhsTTF, RhsTTF>(f, g);
}

template <typename LhsTTF, typename RhsTTF>
TTFWrapper<typename LhsTTF::TTFImpl> merge(const LhsTTF& f, const RhsTTF& g)
{
    typename LhsTTF::UndercutDescriptor dummy_descr;

    return
        merge_impl
            <TTFWrapper<typename LhsTTF::TTFImpl>, LhsTTF, RhsTTF, typename LhsTTF::UndercutDescriptor>
            (f, g, dummy_descr);
}

template <typename LhsTTF, typename RhsTTF, typename UndercutDescr>
TTFWrapper<typename LhsTTF::TTFImpl> merge(const LhsTTF& f, const RhsTTF& g, UndercutDescr& descr)
{
    return merge_impl<TTFWrapper<typename LhsTTF::TTFImpl>, LhsTTF, RhsTTF, UndercutDescr>(f, g, descr);
}


enum class RelativeTtfPosition { UNDEFINED = 0, F_ABOVE = 1, G_ABOVE = 2 };


template <typename LhsTTF, typename RhsTTF>
inline std::deque<std::pair<double,RelativeTtfPosition>>
    analyze_relative_position(const LhsTTF& f, const RhsTTF& g)
{
    assert( LhsTTF::TTFImpl::PERIOD() == RhsTTF::TTFImpl::PERIOD() );
    const double PERIOD = LhsTTF::PERIOD();

    std::deque<std::pair<double,RelativeTtfPosition>> result;

    if ( f.get_max() <= g.get_min() )
    {
         result.push_back( std::make_pair(0.0, RelativeTtfPosition::G_ABOVE) );
         return result;
    }

    if ( f.get_min() >= g.get_max() )
    {
         result.push_back( std::make_pair(0.0, RelativeTtfPosition::F_ABOVE) );
         return result;
    }

    RelativeTtfPosition old_rel_pos = RelativeTtfPosition::UNDEFINED;
    RelativeTtfPosition rel_pos = RelativeTtfPosition::UNDEFINED;

    int i = -1;
    int j = -1;

    if ( eq(f[-1].x, g[-1].x) )
    {
        if ( neq(g[-1].y, f[-1].y) )
        {
            if ( g[-1].y < f[-1].y ) rel_pos = RelativeTtfPosition::F_ABOVE;
            else if ( g[-1].y > f[-1].y ) rel_pos = RelativeTtfPosition::G_ABOVE;
        }
    }
    else if ( g[-1].x < f[-1].x )
    {
        if ( clockwise(g[-1], f[-1], g[0]) ) rel_pos = RelativeTtfPosition::F_ABOVE;
        else if ( counter_clockwise(g[-1], f[-1], g[0]) ) rel_pos = RelativeTtfPosition::G_ABOVE;
    }
    else
    {
        assert( f[-1].x < g[-1].x );
        if ( clockwise(f[-1], g[-1], f[0]) ) rel_pos = RelativeTtfPosition::G_ABOVE;
        else if ( counter_clockwise(f[-1], g[-1], f[0]) ) rel_pos = RelativeTtfPosition::F_ABOVE;
    }

    const int n = std::max(1, (int) f.size());
    const int m = std::max(1, (int) g.size());

    while ( i < n && j < m )
    {
        old_rel_pos = rel_pos;
        if ( eq(f[i+1].x, g[j+1].x) )
        {
            if ( lt(g[j+1].y, f[i+1].y) ) rel_pos = RelativeTtfPosition::F_ABOVE;
            else if ( lt(f[i+1].y, g[j+1].y) ) rel_pos = RelativeTtfPosition::G_ABOVE;
        }
        else if ( f[i+1].x < g[j+1].x )
        {
            if ( clockwise(g[j], f[i+1], g[j+1]) ) rel_pos = RelativeTtfPosition::F_ABOVE;
            if ( counter_clockwise(g[j], f[i+1], g[j+1]) ) rel_pos = RelativeTtfPosition::G_ABOVE;
        }
        else
        {
            assert( g[j+1].x < f[i+1].x );
            if ( clockwise(f[i], g[j+1], f[i+1]) ) rel_pos = RelativeTtfPosition::G_ABOVE;
            if ( counter_clockwise(f[i], g[j+1], f[i+1]) ) rel_pos = RelativeTtfPosition::F_ABOVE;
        }

        if ( old_rel_pos != rel_pos )
        {
            double x;

            if ( eq(f[i].x, g[j].x) && eq(f[i].y, g[j].y) )
                x = f[i].x;
            else if ( colinear(g[j], f[i], g[j+1]) && colinear(f[i], g[j], f[i+1]) )
                x = std::max( f[i].x, g[j].x );
            else if ( colinear(g[j], f[i], g[j+1]) )
                x = f[i].x;
            else if ( colinear(f[i], g[j], f[i+1]) )
                x = g[j].x;
            else
                x = intersection_point(f[i], f[i+1], g[j], g[j+1]).x;
            if ( x >= 0 && x < PERIOD )
                result.push_back( std::make_pair(x, rel_pos) );
        }

          if ( eq(f[i+1].x, g[j+1].x) ) { ++i; ++j; }
          else if ( f[i+1].x < g[j+1].x ) ++i;
          else ++j;
     }

     if ( ! result.empty() && result.front().second == result.back().second )
         result.pop_front();

     if ( result.empty() )
         result.push_back( std::make_pair(0.0, rel_pos) );

    if ( result.size() == 1 && result.front().second == RelativeTtfPosition::UNDEFINED )
         result.front().second = RelativeTtfPosition::G_ABOVE;

    #ifndef NDEBUG
    {
        for ( int i = -1 ; i < (int) result.size() ; ++i )
        {
            assert( KATCH_IMPLIES(result.size() > 1 && i >= 0, result[i].second != RelativeTtfPosition::UNDEFINED) );

            double a;
            double b;
            RelativeTtfPosition pos;

            if ( i < 0 )
            {
                a = result.back().first - PERIOD;
                pos = result.back().second;
            }
            else
            {
                a = result[i].first;
                pos = result[i].second;
            }

            assert( pos == RelativeTtfPosition::F_ABOVE || pos ==  RelativeTtfPosition::G_ABOVE || pos ==  RelativeTtfPosition::UNDEFINED );
            assert( KATCH_IMPLIES(i != 0, pos !=  RelativeTtfPosition::UNDEFINED) );

            if ( i + 1 == result.size() ) b = result[0].first + PERIOD;
            else b = result[i+1].first;

            const int N_TESTS = 100;

            for ( int j = 0 ; j < N_TESTS ; ++j )
            {
                assert(a < b);
                double x = util::random(a, b);

                double fx = f.eval(x);
                double gx = g.eval(x);

                if ( pos == RelativeTtfPosition::F_ABOVE ) assert( gx <= fx + 0.0001 );
                if ( pos == RelativeTtfPosition::G_ABOVE ) assert( fx <= gx + 0.0001 );

                if ( gx < fx - 0.0001 ) assert(  RelativeTtfPosition::F_ABOVE == pos );
                if ( fx < gx - 0.0001 ) assert(  RelativeTtfPosition::G_ABOVE == pos );
            }

        }
    }
    #endif

    return result;
}


}

#endif /* KATCH_TTF_WRAPPER_H_ */
