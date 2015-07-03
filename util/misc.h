/*
 * katch/util/misc.h
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

#ifndef KATCH_MISC_H_
#define KATCH_MISC_H_

#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <chrono>
#include <functional>
#include <limits>
#include <vector>

#include "datastr/graph/basic.h"

#define KATCH_IMPLIES(a, b) ((!(a)) || (b))
#define KATCH_EQUIV(a, b) (KATCH_IMPLIES((a),(b)), KATCH_IMPLIES((b),(a)))

#define KATCH_BE_VERBOSE true

#define KATCH_STATUS(s) if (KATCH_BE_VERBOSE) std::cout << "STATUS [" << __FILE__ << ":" << __LINE__  << "] " << s << std::flush
#define KATCH_WARNING(s) if (KATCH_BE_VERBOSE) std::cerr << "WARNING [" << __FILE__ << ":" << __LINE__  << "] " << s << std::flush
#define KATCH_ERROR(s) if (KATCH_BE_VERBOSE) std::cerr << "ERROR [" << __FILE__ << ":" << __LINE__  << "] " << s << std::flush

#define KATCH_CONTINUE_STATUS(s) if (KATCH_BE_VERBOSE) std::cout << s << std::flush
#define KATCH_CONTINUE_ERROR(s) if (KATCH_BE_VERBOSE) std::cerr << s << std::flush
#define KATCH_CONTINUE_WARNING(s) if (KATCH_BE_VERBOSE) std::cerr << s << std::flush

namespace katch
{
namespace util
{

// Check whether a range (of "less" and "equal_to" comparable elements) contains duplicates
template
<
    typename Iterator,
    typename LtFn = std::less<typename std::iterator_traits<Iterator>::value_type>,
    typename EqFn = std::equal_to<typename std::iterator_traits<Iterator>::value_type>
>
bool has_duplicates
(
        const Iterator& begin,
        const Iterator& end,
        const LtFn& lt = std::less<typename std::iterator_traits<Iterator>::value_type>(),
        const EqFn& eq = std::equal_to<typename std::iterator_traits<Iterator>::value_type>()
)
{
    std::vector<typename std::iterator_traits<Iterator>::value_type> vec(begin, end);

    std::sort(vec.begin(), vec.end(), lt);
    return std::adjacent_find(vec.begin(), vec.end(), eq) != vec.end();
}

// Check whether a vector (of "less" and "equal_to" comparable elements) contains duplicates
template <typename T>
bool has_duplicates(const std::vector<T>& vec)
{
    return has_duplicates(vec.begin(), vec.end());
}


// Remove duplicates from a vector of "less" and "equal_to" comparable elements
template <typename T, typename LtFn = std::less<T>, typename EqFn = std::equal_to<T>>
void remove_duplicates(std::vector<T>& vec, const LtFn& lt = std::less<T>(), const EqFn& eq = std::equal_to<T>())
{
    std::sort(vec.begin(), vec.end(), lt);
    auto new_end = std::unique(vec.begin(), vec.end(), eq);

    vec.resize( std::distance(vec.begin(), new_end) );
    assert( ! has_duplicates(vec) );
}

template <typename Graph, typename Iterator>
std::vector<NodeIterator> neighborhood(const Graph& graph, Iterator begin, Iterator end)
{
    std::vector<NodeIterator> result;

    for ( auto it = begin ; it != end ; ++it )
    {
        const auto& x = *it;
        std::vector<NodeIterator> neighbors_of_x;

        for ( auto e = graph.out_edges_begin(x) ; e != graph.out_edges_end(x) ; ++e )
            neighbors_of_x.push_back(graph.get_target(e));

        for ( auto e = graph.in_edges_begin(x) ; e != graph.in_edges_end(x) ; ++e )
            neighbors_of_x.push_back(graph.get_source(e));

        remove_duplicates(neighbors_of_x);
        result.insert(result.end(), neighbors_of_x.begin(), neighbors_of_x.end());
    }

    return result;
}

// Compute the position of the most significant bit (where the right most position is 0)
inline uint32_t msb(const uint32_t x)
{
    int digit = std::numeric_limits<uint32_t>::digits - 1;
    while( (x & (1 << digit)) == 0 ) --digit;

    assert( digit >= 0 );
    return digit;
}

// Compute x mod m
inline double modulo(const double x, const double m)
{
    if (x >= 0)
    {
        if ( x < m ) return x;
        if ( x < m+m ) return x - m;
    }

    double result = fmod(x, m);
    if (result < 0) result += m;

    assert( 0 <= result );
    assert( result < m );
    assert( fabs(result - (x - m * floor(x/m))) < 0.000001 );

    return result;
}

// compute x mod m
inline int modulo(const int x, const int m)
{
    int result = x % m;
    if (result < 0) result += m;

    assert(result >= 0);
    assert( result == (x - m * floor((double)x/(double)m)) );

    return result;
}


double random(const double a, const double b)
{
    assert( a < b );
    double result =  (double(rand()) / double(RAND_MAX)) * (b - a) + a;

    assert( a <= result + 0.000001 );
    assert( result <= b + 0.000001 );
    return result;
}

using TimePoint = std::chrono::high_resolution_clock::time_point;

inline TimePoint time_stamp()
{
    return std::chrono::high_resolution_clock::now();
}

inline double get_duration_in_seconds(const TimePoint& from, const TimePoint& to)
{
    assert( from <= to );
    return std::chrono::duration_cast<std::chrono::duration<double>>(to - from).count();
}

} /* namespace util */
} /* namespace katch */

#endif /* KATCH_MISC_H_ */
