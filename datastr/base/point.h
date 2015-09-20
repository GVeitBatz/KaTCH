/*
 * katch/datastr/base/point.h
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

#ifndef KATCH_POINT_H_
#define KATCH_POINT_H_

#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>

#include "datastr/base/double.h"

namespace katch
{

struct Point
{
    double x;
    double y;

    Point()
    : x(std::numeric_limits<double>::max()),
      y(std::numeric_limits<double>::max())
    {}

    Point(double xx, double yy) noexcept
    : x(xx), y(yy)
    {}

    friend std::ostream& operator<<(std::ostream& os, const Point& p)
    {
        os << "(" << p.x << ", " << p.y << ")";
        return os;
    }

    friend bool operator== (const Point& lhs, const Point& rhs)
    {
        return lhs.x == rhs.x && lhs.y == rhs.y;
    }

    friend bool operator!= (const Point& lhs, const Point& rhs)
    {
        return ! (lhs == rhs);
    }
};

inline Point operator- (const Point& p, const Point& q) noexcept
{
    return Point(p.x - q.x, p.y - q.y);
}

inline Point operator+ (const Point& p, const Point& q) noexcept
{
    return Point(p.x + q.x, p.y + q.y);
}

inline Point operator* (const double s, const Point& p) noexcept
{
    return Point(s * p.x, s * p.y);
}

inline double perp_dot_product(const Point& p, const Point& q) noexcept
{
    return p.x * q.y - q.x * p.y;
}

//inline int ccw(const Point& a, const Point& b, const Point& c, const double perp_dot_epsilon = 0.0000000001) noexcept
inline int ccw(const Point& a, const Point& b, const Point& c, const double perp_dot_epsilon = 0.000001) noexcept
{
    if ( fabs(a.x - b.x) < 0.000001 && fabs(a.y - b.y) < 0.000001 ) return 0;
    if ( fabs(a.x - c.x) < 0.000001 && fabs(a.y - c.y) < 0.000001 ) return 0;
    if ( fabs(b.x - c.x) < 0.000001 && fabs(b.y - c.y) < 0.000001 ) return 0;

    double x = perp_dot_product(c - a, b - a);

    if ( fabs(x) < perp_dot_epsilon ) return 0;

    if ( x > 0 ) return 1;
    return -1;
}

// checks whether line segments (ab) and (pq) intersect in exactly one point (without border points)
inline bool intersect(const Point& a, const Point& b, const Point& p, const Point& q) noexcept
{
    if ( ccw(a,b,p) == 0 ) return false;
    if ( ccw(a,b,q) == 0 ) return false;
    if ( ccw(p,q,a) == 0 ) return false;
    if ( ccw(p,q,b) == 0 ) return false;

    if ( ccw(a,b,p) == ccw(a,b,q) ) return false;
    if ( ccw(p,q,a) == ccw(p,q,b) ) return false;

    return true;
}

inline Point intersection_point(const Point& a, const Point& b, const Point& p, const Point& q) noexcept
{
    assert( perp_dot_product(p - q, b - a) != 0 );
    return a + (perp_dot_product(p - q, p - a) / perp_dot_product(p - q, b - a)) * (b - a);
}

// computes intersection point of line segment (ab) with a.x < b.x and the horizontal line given by the equation y = c
inline Point intersection_point(const Point& a, const Point& b, const double c)
{
    assert( std::min(a.y, b.y) <= c + 0.000001 );
    assert( c <= std::max(a.y, b.y) + 0.000001 );

    const Point p = Point(a.x, c);
    const Point q = Point(b.x, c);

    if ( fabs(p.x - q.x) <= 0.000001 )
    {
        assert( std::min(p.y, q.y) <= c + 0.000001 );
        assert( c <= std::max(p.y, q.y) + 0.000001 );

        return Point((p.x + q.x) / 2.0, c);
    }

    assert( intersect(a, b, p, q) );
    Point result = intersection_point(a, b, p, q);

    assert( std::min(a.y, b.y) <= result.y + 0.000001 );
    assert( result.y <= std::max(a.y, b.y) + 0.000001 );
    return result;
}

// tells whether (p,q,r) are oriented clockwise in the plain
static inline bool eq(const Point& p, const Point& q) { return eq(p.x, q.x) && eq(p.y, q.y); }

// tells whether (p,q,r) are oriented clockwise in the plain
inline bool clockwise(const Point& p, const Point& q, const Point& r) { return ccw(p,q,r) == 1; }

// tells whether (p,q,r) are colinear in the plain
inline bool colinear(const Point& p, const Point& q, const Point& r) { return ccw(p,q,r) == 0; }

// tells whether (p,q,r) are oriented counterclockwise in the plain
inline bool counter_clockwise(const Point& p, const Point& q, const Point& r) { return ccw(p,q,r) == -1; }

// tells whether the counterclockwise angle spanned by (p,q,r), with apex q, is less than Pi
inline bool upper_angle_lt_pi(const Point& p, const Point& q, const Point& r)
{
    return counter_clockwise(q, p, r);
}

// tells whether the counterclockwise angle spanned by (p,q,r), with apex q, is greater than Pi
inline bool upper_angle_gt_pi(const Point& p, const Point& q, const Point& r)
{
    return clockwise(q, p, r);
}

// tells whether the clockwise angle spanned by (p,q,r), with apex q, is less than Pi
inline bool lower_angle_lt_pi(const Point& p, const Point& q, const Point& r)
{
    return clockwise(q, p, r);
}

// tells whether the clockwise angle spanned by (p,q,r), with apex q, is greater than Pi
inline bool lower_angle_gt_pi(const Point& p, const Point& q, const Point& r)
{
    return counter_clockwise(q, p, r);
}

}


#endif /* KATCH_POINT_H_ */
