/*
 * katch/datastr/base/interval.h
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

#ifndef KATCH_INTERVAL_H_
#define KATCH_INTERVAL_H_

#include <assert.h>
#include <limits>

namespace katch
{

class Interval
{
private:

    double _lower;
    double _upper;

public:

    static const Interval INFTY;

    constexpr Interval()
    : _lower(std::numeric_limits<double>::max()), _upper(std::numeric_limits<double>::max())
    {}

    constexpr Interval(const double lower, const double upper)
    : _lower(lower), _upper(upper)
    {}

    double get_lower() const
    {
        return _lower;
    }

    double get_upper() const
    {
        return _upper;
    }

    friend Interval merge(const Interval& lhs, const Interval& rhs)
    {
        return
            Interval(
                std::min(lhs.get_lower(), rhs.get_lower()),
                std::min(lhs.get_upper(), rhs.get_upper())
            );
    }

    friend bool operator== (const Interval& lhs, const Interval& rhs)
    {
        return lhs._lower == rhs._lower && lhs._upper == rhs._upper;
    }

    friend bool operator!= (const Interval& lhs, const Interval& rhs)
    {
        return ! ( lhs == rhs );
    }

    friend std::ostream& operator<<(std::ostream& os, const Interval& interval)
    {
        os << "[" << interval._lower << "," << interval._upper << "]";
        return os;
    }
};

const Interval Interval::INFTY(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());

}


#endif /* 'KATCH_INTERVAL_H_ */
