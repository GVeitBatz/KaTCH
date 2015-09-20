/*
 * katch/datastr/base/double.h
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

#ifndef KATCH_DOUBLE_H_
#define KATCH_DOUBLE_H_

namespace katch
{

constexpr double EPSILON = 0.000001;

inline bool le(const double& x, const double& y) { return x <= y + EPSILON; }
inline bool lt(const double& x, const double& y) { return x + EPSILON < y; }
inline bool eq(const double& x, const double& y) { return fabs(x - y) <= EPSILON; }
inline bool neq(const double& x, const double& y) { return ! eq(x, y); }
inline bool gt(const double& x, const double& y) { return lt(y, x); }
inline bool ge(const double& x, const double& y) { return le(y, x); }

}

#endif /* KATCH_DOUBLE_H_ */
