/*
 * katch/datastr/base/pwl_ttf.h
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

#ifndef KATCH_PWL_TTF_H_
#define KATCH_PWL_TTF_H_

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <limits>
#include <vector>
#include <iomanip>

#include "datastr/base/double.h"
#include "datastr/base/point.h"
#include "util/misc.h"

namespace katch
{

template <int Period, int BucketSize = 8>
class PwlTTF
{

public:

    static constexpr double PERIOD() { return Period; }
    static constexpr int BUCKET_SIZE() { return BucketSize; }

    struct UndercutDescriptor
    {
        bool _f_undercuts_strictly;
        bool _g_undercuts_strictly;
    };

private:

    static constexpr int N_TEST_CASES = 100;

    static inline void append_point(std::vector<Point>& points, const Point& p)
    {
        assert( p.y >= 0 );
        assert(KATCH_IMPLIES( ! points.empty(), le(points.back().x, p.x) ));

        if ( ! points.empty() && eq(points.back().x, p.x) )
        {
            if ( eq(points.back().y, p.y) ) return;

            KATCH_WARNING("Detected discontinuity " << points.back() << ", " << p);

            Point pp( std::max(points.back().x + EPSILON, p.x + EPSILON), p.y );
            points.back().x = std::min( points.back().x, p.x );

            assert( points.back().x < pp.x );
            points.push_back(pp);

            return;
        }

        if ( ! points.empty() && p.x < points.back().x )
        {
            KATCH_WARNING("Appending  point " << p << " left of " << points.back());

            Point pp(points.back().x, p.y);
            points.back().x = p.x;

            assert( points.back().x < pp.x );
            points.push_back(pp);

            return;
        }

        if ( points.size() > 1 )
        {
            if ( colinear( *(points.end() - 2), *(points.end() - 1), p ) )
            {
                points.back() = p;
                return;
            }
        }

        assert(KATCH_IMPLIES( ! points.empty(), points.back().x < p.x ));
        points.push_back(p);
    }

    template <typename PointIterator>
    static size_t compute_bucket_shift(const PointIterator begin, const PointIterator end)
    {
        uint32_t n_points =  begin < end ? end - begin : 1;
        return util::msb( (uint32_t(Period) * uint32_t(BucketSize)) / n_points );
    }

    static uint32_t compute_n_buckets(const uint32_t bucket_shift)
    {
        return std::max(uint32_t(1), uint32_t(Period) >> bucket_shift);
    }

    // yields the number of the bucket where x belongs to
    inline int get_bucket(const double& x) const
    {
        assert( 0 <= x && x <= PERIOD() );

        int bucket( (static_cast<size_t>(x) >> _bucket_shift) );
        if ( bucket > 0 ) --bucket;

        assert( bucket < _buckets.size() );
        assert( bucket >= 0 );

        return bucket;
    }

    void set_up_buckets_and_min_max()
    {
        int bucket = -1;
        for ( auto it = _points.begin() ; it != _points.end() ; ++it )
        {
            while ( get_bucket(it->x) > bucket )
            {
                assert( bucket+1 < _buckets.size() );
                _buckets[++bucket] = it - _points.begin();
            }

            _min = std::min(_min, it->y);
            _max = std::max(_max, it->y);

        }

        while ( ++bucket < _buckets.size() ) _buckets[bucket] = _points.size() - 1;
    }


    int compute_first_index_with_greater_or_equal(const double& x) const
    {
        double xx = util::modulo(x, PERIOD());
        assert( 0 <= xx && xx < PERIOD() );
        assert( _buckets[0] == 0 );

        const int offset = _points.size() * static_cast<size_t>( floor(x / PERIOD()) );

        if ( xx <= _points[0].x )
        {
            assert( (*this)[offset-1].x < x );
            assert( (*this)[offset].x >= x );
            return offset;
        }

        const int b = get_bucket(xx);
        assert( 0 <= b && b < _buckets.size() );
        assert( KATCH_IMPLIES( _buckets[b] > 0, le(_points[_buckets[b] - 1].x, xx) ) );

        int index = _buckets[b];
        assert( 0 <= index && index < _points.size() );
        assert( 0 <= _points[index].x && _points[index].x < PERIOD() );

        while ( _points[index].x < xx )
        {
            ++index;

            if ( index >= _points.size() )
            {
                assert( (*this)[_points.size() + offset - 1].x < x );
                assert( (*this)[_points.size() + offset].x >= x );
                return _points.size() + offset;
            }

            assert( 0 <= index && index < _points.size() );
            assert( 0 <= _points[index].x && _points[index].x < PERIOD() );
        }

        assert( index < _points.size() );
        assert( (*this)[index - 1 + offset].x < x );
        assert( (*this)[index + offset].x >= x );

        return index + offset;
    }

    bool dbg_check_points_valid()
    {
        if ( _points.empty() ) return false;

        for ( auto it = _points.begin() + 1 ; it != _points.end() ; ++it )
        {
            if ( (it-1)->x >= it->x ) return false;
            if ( (it-1)->y < 0 ) return false;

            if ( eq((it-1)->x, it->x) ) KATCH_WARNING("Detected nearby points");
        }

        return true;
    }

    bool dbg_check_eval(const double time, const double result) const
    {
        const double x = util::modulo(time, PERIOD());

        Point p, q;

        if ( lt(x, _points.front().x) )
        {
            p = _points.back() + Point(-PERIOD(), 0);
            q = _points.front();
        }
        else if ( gt(x, _points.back().x) )
        {
            p = _points.back();
            q = _points.front() + Point(PERIOD(), 0);
        }
        else
        {
            for ( size_t i = 0 ; i + 1 < _points.size() ; ++i )
            {
                if ( eq(_points[i].x, x) ) return _points[i].y;

                if ( gt( _points[i+1].x, x) )
                {
                    assert( le(_points[i].x, x) );
                    p = _points[i];
                    q = _points[i+1];

                    break;
                }
            }
        }

        double m = (q.y - p.y) / (q.x - p.x);
        return eq(m * (x - p.x) + p.y, result);
    }

    // check whether h = g*f
    static bool dbg_check_link(PwlTTF const * const g_ptr, PwlTTF const * const f_ptr, PwlTTF const * const h_ptr)
    {
        const PwlTTF& f = *f_ptr;
        const PwlTTF& g = *g_ptr;
        const PwlTTF& h = *h_ptr;

        std::vector<double> x_values;

        for ( int i = 0 ; i < g.get_n_points() ; ++i ) x_values.push_back(g[i].x);
        for ( int i = 0 ; i < f.get_n_points() ; ++i ) x_values.push_back(f[i].x);
        for ( int i = 0 ; i < h.get_n_points() ; ++i ) x_values.push_back(h[i].x);
        for ( int i = 0 ; i < N_TEST_CASES ; ++i ) x_values.push_back( util::random(-2.0*PERIOD(), 2.0*PERIOD()) );

        x_values.push_back(0.0);
        x_values.push_back(PERIOD());
        x_values.push_back(-PERIOD());
        x_values.push_back(2.0*PERIOD());
        x_values.push_back(-2.0*PERIOD());

        for ( auto it = x_values.begin() ; it != x_values.end() ; ++it )
        {
            const double x = *it;
            const double result = h.eval(x);
            const double reference_result = g.eval(f.eval(x) + x) + f.eval(x);

            if ( fabs(result - reference_result) >= 0.0001 ) return false;
        }

        return true;
    }

    // check whether h = min(f,g)
    static bool dbg_check_min(PwlTTF const * const f_ptr, PwlTTF const * const g_ptr, PwlTTF const * const h_ptr)
    {
        double DBG_EPSILON = 0.0001;

        const PwlTTF& f = *f_ptr;
        const PwlTTF& g = *g_ptr;
        const PwlTTF& h = *h_ptr;

        std::vector<double> x_values;

        x_values.push_back(0.0);
        x_values.push_back(PERIOD());

        for ( int i = 0 ; i < g.get_n_points() ; ++i ) x_values.push_back(g[i].x);
        for ( int i = 0 ; i < f.get_n_points() ; ++i ) x_values.push_back(f[i].x);
        for ( int i = 0 ; i < h.get_n_points() ; ++i ) x_values.push_back(h[i].x);
        for ( int i = 0 ; i < N_TEST_CASES ; ++i ) x_values.push_back( util::random(-2.0*PERIOD(), 2.0*PERIOD()) );

        x_values.push_back(0.0);
        x_values.push_back(PERIOD());
        x_values.push_back(-PERIOD());
        x_values.push_back(2.0*PERIOD());
        x_values.push_back(-2.0*PERIOD());

        for ( auto it = x_values.begin() ; it != x_values.end() ; ++it )
        {
            double x = *it;
            double result = h.eval(x);
            double reference_result = std::min(f.eval(x), g.eval(x));

            if ( fabs(result - reference_result) > DBG_EPSILON )
                return false;
        }

        return true;
    }

    // check whether UndercutDescriptor setup by merging of TTFs h = merge(f,g) is correct
    static bool dbg_check_undercut_descr(PwlTTF const * const f_ptr, PwlTTF const * const g_ptr, PwlTTF const * const h_ptr, const UndercutDescriptor& descr)
    {
        double DBG_EPSILON = 0.0001;

        const PwlTTF& f = *f_ptr;
        const PwlTTF& g = *g_ptr;
        const PwlTTF& h = *h_ptr;

        bool f_undercuts_g_strictly = false;
        bool g_undercuts_f_strictly = false;

        for ( int i = 0 ; i < f.get_n_points() ; ++i )
        {
            if ( f[i].y + DBG_EPSILON < g.eval(f[i].x) ) f_undercuts_g_strictly = true;
            if ( g.eval(f[i].x) + DBG_EPSILON < f[i].y ) g_undercuts_f_strictly = true;
        }

        for ( int i = 0 ; i < g.get_n_points() ; ++i )
        {
            if ( g[i].y + DBG_EPSILON < f.eval(g[i].x ) ) g_undercuts_f_strictly = true;
            if ( f.eval(g[i].x) + DBG_EPSILON < g[i].y ) f_undercuts_g_strictly = true;
        }

        if ( ! KATCH_IMPLIES(f_undercuts_g_strictly, descr._f_undercuts_strictly) )
        {
            KATCH_WARNING( \
                    "Unclear situation: f_undercuts_g_strictly = " << f_undercuts_g_strictly \
                    << ", descr._f_undercuts_strictly = " << descr._f_undercuts_strictly << "\n");

            return false;
        }

        if ( ! KATCH_IMPLIES(g_undercuts_f_strictly, descr._g_undercuts_strictly) )
        {
            KATCH_WARNING( \
                    "Unclear situation: g_undercuts_f_strictly = " << g_undercuts_f_strictly \
                    << ", descr._g_undercuts_strictly = " << descr._g_undercuts_strictly << "\n");

            return false;
        }

        std::vector<double> x_values;

        for ( int i = 0 ; i < g.get_n_points() ; ++i ) x_values.push_back(g[i].x);
        for ( int i = 0 ; i < f.get_n_points() ; ++i ) x_values.push_back(f[i].x);
        for ( int i = 0 ; i < h.get_n_points() ; ++i ) x_values.push_back(h[i].x);
        for ( int i = 0 ; i < N_TEST_CASES ; ++i ) x_values.push_back( util::random(-2.0*PERIOD(), 2.0*PERIOD()) );

        x_values.push_back(0.0);
        x_values.push_back(PERIOD());
        x_values.push_back(-PERIOD());
        x_values.push_back(2.0*PERIOD());
        x_values.push_back(-2.0*PERIOD());

        for ( auto it = x_values.begin() ; it != x_values.end() ; ++it )
        {
            double x = *it;

            if ( ! KATCH_IMPLIES( f.eval(x) + DBG_EPSILON < g.eval(x), descr._f_undercuts_strictly ) ) return false;
            if ( ! KATCH_IMPLIES( g.eval(x) + DBG_EPSILON < f.eval(x), descr._g_undercuts_strictly ) ) return false;
        }

        return true;
    }

    const std::vector<Point>& get_ttf_points() const
    {
        return _points;
    }

    uint32_t _bucket_shift;
    std::vector<Point> _points;
    std::vector<int> _buckets;
    double _min;
    double _max;

    PwlTTF(std::vector<Point>&& points) noexcept
    : _bucket_shift( compute_bucket_shift(points.begin(), points.end()) ),
      _points( std::move(points) ),
      _buckets( compute_n_buckets(_bucket_shift), 0 ),
      _min( std::numeric_limits<double>::max() ),
      _max( std::numeric_limits<double>::min() )
    {
        set_up_buckets_and_min_max();

        assert( dbg_check_points_valid() );
    }

public:

    PwlTTF() noexcept
    : _bucket_shift(0), _points(), _buckets(),
      _min(std::numeric_limits<double>::max()),
      _max(std::numeric_limits<double>::min())
    {}

    template <typename Iterator>
    PwlTTF(Iterator points_begin, Iterator points_end) noexcept
    : _bucket_shift(0),
      _points(),
      _buckets(),
      _min( std::numeric_limits<double>::max() ),
      _max( std::numeric_limits<double>::min() )
    {
        for ( auto it = points_begin ; it != points_end ; ++it )
            append_point( _points, *it );

        _bucket_shift = compute_bucket_shift(_points.begin(), _points.end());
        _buckets.assign(compute_n_buckets(_bucket_shift), 0);

        set_up_buckets_and_min_max();

        assert( dbg_check_points_valid() );
    }


    PwlTTF(const double c) noexcept
    : PwlTTF( std::vector<Point>(1, Point(0, c)) )
    {}

    bool is_const() const noexcept
    {
        return _points.size() == 1;
    }

    bool is_fifo() const noexcept
    {
        for ( int i = 0 ; i < _points.size() - 1 ; i++)
        {
            Point p( _points[i] );
            Point q( _points[i+1] );

            if ( lt(q.x + q.y, p.x + p.y) ) return false;
        }

        Point p( _points.back() );
        Point q( Point(_points.front().x + PERIOD(), _points.front().y) );

        if ( lt(q.x + q.y, p.x + p.y) ) return false;

        return true;
    }


    // get i-th point of TTF with i from Z
    Point operator [] (const int i) const noexcept
    {
        assert( ! _points.empty() );

        if ( i >= 0 && i < _points.size() )
            return _points[i];

        if ( i >= -_points.size() && i < 0 )
            return _points[i + _points.size()] + Point(-PERIOD(), 0);

        if ( i >= _points.size() && i < _points.size() + _points.size() )
            return _points[i - _points.size()] + Point(PERIOD(), 0);

        int i_mod = i % _points.size();
        if ( i_mod < 0 ) i_mod += _points.size();

        const double offset = PERIOD() * double(i / _points.size());

        assert( 0 <= i_mod );
        assert( i_mod < _points.size() );
        return _points[i_mod] + Point(offset, 0);
    }

    double get_const_value() const noexcept
    {
        assert( is_const() );
        return _points.back().y;
    }

    size_t get_n_points() const noexcept
    {
        return _points.size();
    }

    size_t size() const noexcept
    {
        return get_n_points();
    }

    double get_min() const noexcept
    {
        assert( this );

        assert( _min ==
            std::min_element(_points.begin(), _points.end(),
                                [](const Point& p, const Point& q) { return p.y < q.y; })->y
        );

        return _min;
    }

    double get_max() const noexcept
    {
        assert( _max ==
            std::max_element(_points.begin(), _points.end(),
                                [](const Point& p, const Point& q) { return p.y < q.y; })->y
        );

        return _max;
    }

    // evaluate the TTF for a given departure time; that is, compute f(time)
    double eval(const double time) const noexcept
    {
        if ( is_const() ) return _points.back().y;

        const double x = util::modulo(time, PERIOD());
        assert( x >= 0 );
        assert( x <= PERIOD() );

        const int bucket = get_bucket(x);
        assert( 0 <= bucket );
        assert( bucket < _buckets.size() );

        assert( _buckets[bucket] >= 0 );
        assert( _buckets[bucket] < _points.size() );
        assert(KATCH_IMPLIES( _buckets[bucket] > 0, le(_points[_buckets[bucket] - 1].x, x) ));

        for ( int i = _buckets[bucket] ; i != _points.size() ; ++i )
        {
            assert( 0 <= i && i < _points.size() );

            if ( le(x, _points[i].x) )
            {
                 if ( eq(x, _points[i].x) ) return _points[i].y;

                 Point p( i > 0 ? _points[i-1] : Point(_points.back().x - PERIOD(), _points.back().y) );
                 Point q( _points[i] );

                 assert( le(p.x, x) );
                 assert( le(x, q.x) );

                 double m = (q.y - p.y) / (q.x - p.x);
                 double result = m * (x - p.x) + p.y;

                return result;
            }
        }

        assert( le(_points.back().x, x) );
        assert( le(x, _points.front().x + PERIOD()) );

        Point q( Point(_points.front().x + PERIOD(), _points.front().y) );
        Point p( _points.back() );

        double m = (q.y - p.y) / (q.x - p.x);
        double result = m * (x - p.x) + p.y;

        return result;
    }

    // compute f+c, where f is a TTF and c is a constant
    PwlTTF* add(const double c) const noexcept
    {
        const PwlTTF& f = *this;
        std::vector<Point> result_points;
        result_points.reserve(f.size() + 2);

        for ( int i = 0 ; i < f.get_n_points() ; ++i )
            append_point( result_points, f[i] + Point(0, c) );

        PwlTTF* result_ptr = new PwlTTF( std::move(result_points) );
        return result_ptr;
    }

    // compute s*f, where f is a TTF and s is a non-negative real number
    PwlTTF* scale(const double s) const noexcept
    {
        assert( s >= 0 );

        const PwlTTF& f = *this;
        std::vector<Point> result_points;
        result_points.reserve(f.size() + 2);

        for ( int i = 0 ; i < f.get_n_points() ; ++i )
            append_point( result_points, Point(f[i].x, s * f[i].y) );

        PwlTTF* result_ptr = new PwlTTF( std::move(result_points) );
        return result_ptr;
    }

    // computes the TTF g*f of the linked path u -f-> v -g-> w; that is, g "after" f
    friend inline PwlTTF* link(PwlTTF const * const g_ptr, PwlTTF const * const f_ptr) noexcept
    {
        assert( g_ptr );
        assert( f_ptr );
        assert( f_ptr->is_fifo() );

        if ( g_ptr->is_const() && f_ptr->is_const() )
            return new PwlTTF( g_ptr->get_const_value() + f_ptr->get_const_value() );

        if ( f_ptr->is_const() )
            return link(g_ptr, f_ptr->get_const_value());

        if ( g_ptr->is_const() )
            return link(g_ptr->get_const_value(), f_ptr);

        const PwlTTF& f( *f_ptr );
        const PwlTTF& g( *g_ptr );

        std::vector<Point> result_points;
        result_points.reserve(f.size() + g.size() + 10);

        int i = g.compute_first_index_with_greater_or_equal(f.eval(0));
        int j = 0;

        while( true )
        {
            double x, y;

            if ( eq( g[i].x, f[j].x + f[j].y ) )
            {
                x = f[j].x;
                y = g[i].y + f[j].y;

                ++i; ++j;
            }
            else if ( g[i].x < f[j].x + f[j].y )
            {
                assert( lt(g[i].x, f[j].x + f[j].y) );

                double m_arr_f_inverse = (f[j].x - f[j-1].x) / (f[j].x + f[j].y - f[j-1].x - f[j-1].y);
                x = m_arr_f_inverse * (g[i].x - f[j-1].x - f[j-1].y) + f[j-1].x;
                y = g[i].x + g[i].y - x;

                ++i;
            }
            else
            {
                assert( lt(f[j].x + f[j].y, g[i].x) );

                x = f[j].x;
                double m_g = (g[i].y - g[i-1].y) / (g[i].x - g[i-1].x);
                y = g[i-1].y + m_g * (f[j].x + f[j].y - g[i-1].x) + f[j].y;

                ++j;
            }

            if ( lt(PERIOD(), x) ) break;

            x = std::min(x, PERIOD());
            x = std::max(x, 0.0);

            append_point(result_points, Point(x,y));
        }

        if ( result_points.size() > 1 )
        {
            size_t last = result_points.size() - 1;

            if ( colinear(result_points[last-1], result_points[last], result_points[0] + Point(PERIOD(), 0)) )
                result_points.pop_back();
        }

        result_points.shrink_to_fit();
        PwlTTF* result_ptr = new PwlTTF( std::move(result_points) );

        assert( dbg_check_link(g_ptr, f_ptr, result_ptr) );
        return result_ptr;
    }

    // computes the TTF f*c of the linked path u -f-> v -c-> w; that is, f "after" c where c is a constant function
    friend inline PwlTTF* link(PwlTTF const * const f_ptr, const double c) noexcept
    {
        const PwlTTF& f = *f_ptr;
        const double c_mod = util::modulo(c, PERIOD());

        int b = f.get_bucket(c_mod);
        int k = f._buckets[b];

        while ( k < f.get_n_points() && lt(f[k].x, c_mod) ) ++k;

        std::vector<Point> result_points;
        result_points.reserve(f.size() + 2);

        for ( int i = k ; i < f.get_n_points() ; ++i )
            append_point( result_points, f[i] + Point(-c_mod, c) );

        for ( int i = 0 ; i < k ; ++i )
            append_point( result_points, f[i] + Point(PERIOD() - c_mod, c) );

        assert( ! result_points.empty() );
        result_points.shrink_to_fit();
        PwlTTF* result_ptr = new PwlTTF( std::move(result_points) );

        PwlTTF constant_ttf(c);
        assert( dbg_check_link(f_ptr, &constant_ttf, result_ptr) );
        return result_ptr;
    }

    // computes the TTF c*f = f+c of the linked path u -f-> v -c-> w; that is, c "after" f where c is a constant function
    friend inline PwlTTF* link(const double c, PwlTTF const * const f_ptr) noexcept
    {
        PwlTTF* result_ptr = f_ptr->add(c);

        PwlTTF constant_ttf(c);
        assert( dbg_check_link(&constant_ttf, f_ptr, result_ptr) );
        return result_ptr;
    }

    friend inline PwlTTF* merge(PwlTTF const * const f_ptr, PwlTTF const * const g_ptr, UndercutDescriptor& descr) noexcept
    {
        descr = { false, false };

        const PwlTTF& f = *f_ptr;
        const PwlTTF& g = *g_ptr;

        if ( lt(f.get_max(), g.get_min()) )
        {
            descr._f_undercuts_strictly = true;
            return new PwlTTF(f);
        }
        else if ( lt(g.get_max(), f.get_min()) )
        {
            descr._g_undercuts_strictly = true;
            return new PwlTTF(g);
        }

        std::vector<Point> result_points;
        result_points.reserve(2*f.size() + 2*g.size() + 2);

        int i = 0, j = 0;

        while ( i < f.get_n_points() || j < g.get_n_points() )
        {
            if ( intersect(f[i-1], f[i], g[j-1], g[j]) )
            {
                const Point s = intersection_point(f[i-1], f[i], g[j-1], g[j]);
                if (s.x >= 0) append_point(result_points, s);

                descr._f_undercuts_strictly = true;
                descr._g_undercuts_strictly = true;
            }

            if ( eq(f[i].x, g[j].x) )
            {
                if ( eq(f[i].y, g[j].y) )
                {
                    append_point(result_points, f[i]);

                    if ( counter_clockwise(g[j-1], f[i-1], f[i]) || counter_clockwise(f[i], f[i+1], g[j+1]) )
                        descr._f_undercuts_strictly = true;

                    if ( clockwise(g[j-1], f[i-1], f[i]) || clockwise(f[i], f[i+1], g[j+1]) )
                        descr._g_undercuts_strictly = true;
                }
                else if ( f[i].y < g[j].y )
                {
                    append_point(result_points, f[i]);

                    assert( lt(f[i].y, g[j].y) );
                    descr._f_undercuts_strictly = true;
                }
                else
                {
                    assert( g[j].y < f[i].y );
                    append_point(result_points, g[j]);
                    descr._g_undercuts_strictly = true;
                }

                i++; j++;
            }
            else if (f[i].x < g[j].x)
            {
                assert( lt(f[i].x, g[j].x) );

                if ( counter_clockwise(g[j-1], f[i], g[j]) )
                {
                    descr._f_undercuts_strictly = true;
                    append_point(result_points, f[i]);
                }
                else if ( colinear(g[j-1], f[i], g[j]) )
                {
                    if ( counter_clockwise(g[j-1], f[i-1], f[i]) || counter_clockwise(f[i], f[i+1], g[j]) )
                    {
                        descr._f_undercuts_strictly = true;
                        append_point(result_points, f[i]);
                    }
                    else if ( result_points.empty() )
                        append_point(result_points, f[i]);

                    if ( clockwise(g[j-1], f[i-1], f[i]) || clockwise(f[i], f[i+1], g[j]) )
                        descr._g_undercuts_strictly = true;
                }

                ++i;
            }
            else
            {
                assert( g[j].x < f[i].x );

                if ( counter_clockwise(f[i-1], g[j], f[i]) )
                {
                    descr._g_undercuts_strictly = true;
                    append_point(result_points, g[j]);
                }
                else if ( colinear(f[i-1], g[j], f[i]) )
                {
                    if ( counter_clockwise(g[j-1], g[j], f[i-1]) || counter_clockwise(g[j], g[j+1], f[i]) )
                    {
                        descr._g_undercuts_strictly = true;
                        append_point(result_points, g[j]);
                    }
                    if ( result_points.empty() )
                        append_point(result_points, g[j]);

                    if ( clockwise(g[j-1], g[j], f[i-1]) || clockwise(g[j], g[j+1], f[i]) )
                        descr._f_undercuts_strictly = true;
                }

                ++j;
            }
        }

        const int n = f.get_n_points();
        const int m = g.get_n_points();

        if ( intersect(f[n-1], f[n], g[m-1], g[m]) )
        {
            const Point s = intersection_point(f[n-1], f[n], g[m-1], g[m]);
            if (s.x < PERIOD()) append_point(result_points, s);

            descr._f_undercuts_strictly = true;
            descr._g_undercuts_strictly = true;
        }

        assert( ! result_points.empty() );
        result_points.shrink_to_fit();
        PwlTTF* result_ptr = new PwlTTF( std::move(result_points) );

        assert( dbg_check_min(f_ptr, g_ptr, result_ptr) );
        assert( dbg_check_undercut_descr(f_ptr, g_ptr, result_ptr, descr) );
        return result_ptr;
    }

    friend inline PwlTTF* merge(PwlTTF const * const f_ptr, const double c, UndercutDescriptor& descr) noexcept
    {
        descr = { false, false };

        const PwlTTF& f = *f_ptr;
        std::vector<Point> result_points;
        result_points.reserve(f.size() + 2);

        if ( lt(c, f.get_max()) ) descr._g_undercuts_strictly = true;
        if ( lt(f.get_min(), c) ) descr._f_undercuts_strictly = true;

        if ( le(c, f.get_min()) ) return new PwlTTF(c);

        for (int i = 0 ; i < f.get_n_points() ; ++i)
        {
            if ( eq(f[i].y, c) )
            {
                if ( lt(f[i-1].y, c) || lt(f[i+1].y, c) )
                {
                    append_point(result_points, Point(f[i].x, c));
                    descr._f_undercuts_strictly = true;
                }
                else if ( result_points.empty() )
                    append_point(result_points, Point(f[i].x, std::min(f[i].y,c)));

                if ( lt(c, f[i-1].y) || lt(c, f[i+1].y) )
                    descr._g_undercuts_strictly = true;
            }
            else if ( f[i].y < c )
            {
                assert( descr._f_undercuts_strictly );

                if ( lt(c, f[i-1].y) )
                {
                    assert( descr._g_undercuts_strictly );
                    const Point s = intersection_point(f[i-1], f[i], c);

                    if ( s.x >= 0 )
                    {
                        assert( le(s.x, PERIOD()) );
                        assert( le(0, s.y) );
                        assert( eq(s.y, c) );

                        append_point(result_points, s);
                    }
                }

                append_point(result_points, f[i]);
            }
            else if ( lt(f[i-1].y, c) )
            {
                assert( f[i].y > c );
                assert( descr._f_undercuts_strictly );
                assert( descr._g_undercuts_strictly );

                const Point s = intersection_point(f[i-1], f[i], c);

                if ( s.x >= 0 )
                {
                    assert( le(s.x, PERIOD()) );
                    assert( le(0, s.y) );
                    assert( eq(s.y, c) );

                    append_point(result_points, s);
                }
            }
        }

        const int n = f.get_n_points();

        if ( (lt( f[n-1].y, c) && lt(c, f[n].y)) || (lt(c, f[n-1].y) && lt(f[n].y, c)) )
        {
            assert( descr._f_undercuts_strictly );
            assert( descr._g_undercuts_strictly );
            const Point s = intersection_point(f[n-1], f[n], c);

            if ( le(s.x, PERIOD()) )
            {
                assert( le(0, s.x) );
                append_point(result_points, s);
            }
        }

        if ( result_points.empty() )
            result_points.push_back(Point(0.0, c));

        result_points.shrink_to_fit();
        PwlTTF* result_ptr = new PwlTTF( std::move(result_points) );

        PwlTTF constant_ttf(c);

        assert( dbg_check_min(f_ptr, &constant_ttf, result_ptr) );
        assert( dbg_check_undercut_descr(f_ptr, &constant_ttf, result_ptr, descr) );
        return result_ptr;
    }

    friend inline PwlTTF* merge(const double c, PwlTTF const * const g_ptr, UndercutDescriptor& descr) noexcept
    {
        UndercutDescriptor descr_swapped;
        PwlTTF* result_ptr = merge(g_ptr, c, descr_swapped);

        descr._g_undercuts_strictly = descr_swapped._f_undercuts_strictly;
        descr._f_undercuts_strictly = descr_swapped._g_undercuts_strictly;

        return result_ptr;
    }

    friend inline PwlTTF* merge(PwlTTF const * const f_ptr, PwlTTF const * const g_ptr) noexcept
    {
        UndercutDescriptor dummy_descr;
        return merge(f_ptr, g_ptr, dummy_descr);
    }

    friend inline PwlTTF* merge(PwlTTF const * const f_ptr, const double c) noexcept
    {
        UndercutDescriptor dummy_descr;
        return merge(f_ptr, c, dummy_descr);
    }

    friend inline PwlTTF* merge(const double c, PwlTTF const * const g_ptr) noexcept
    {
        UndercutDescriptor dummy_descr;
        return merge(c, g_ptr, dummy_descr);
    }

    friend std::ostream& operator<<(std::ostream& os, const PwlTTF& f)
    {
        os << "{ ";

        for ( int i = 0 ; i < f.size() ; ++i )
            os << "(" << f[i].x << "," << f[i].y << ") ";

        os << "}";

        return os;
    }
};

}

#endif /* KATCH_PWL_TTF_H_ */
