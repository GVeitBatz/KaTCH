/*
 * katch/io/edge_info.h
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

#ifndef KATCH_EDGE_INFO_H_
#define KATCH_EDGE_INFO_H_

#include <cassert>
#include <limits>
#include <vector>

#include "datastr/graph/basic.h"

namespace katch
{

template <typename TTFType>
class EdgeInfo
{

public:

    using TTF = TTFType;

private:

    size_t _source;
    size_t _target;
    bool _forward : 1;
    bool _backward : 1;
    TTF _ttf;
    std::vector<MiddleNodeDescr> _middle_node_data;

public:

    EdgeInfo()
    : _source(std::numeric_limits<size_t>::max()),
      _target(std::numeric_limits<size_t>::max()),
      _forward(false),
      _backward(false),
      _ttf(),
      _middle_node_data()
    {}

    size_t get_source() const { return _source; }
    void set_source(const size_t& source) { _source = source; }

    size_t get_target() const { return _target; }
    void set_target(const size_t& target) { _target = target; }

    bool get_forward() const { return _forward; }
    void set_forward(const bool& forward) { _forward = forward; }

    bool get_backward() const { return _backward; }
    void set_backward(const bool& backward) { _backward = backward; }

    void set_ttf(TTF&& ttf) { _ttf = std::move(ttf); }
    void set_ttf(const double c) { _ttf = TTF(c); }
    const TTF& get_ttf() const { return _ttf; }
    TTF& get_ttf() { return _ttf; }

    void set_middle_node_data(std::vector<MiddleNodeDescr>&& middle_node_data)
    {
        _middle_node_data = std::move(middle_node_data);
    }

    const std::vector<MiddleNodeDescr>& get_middle_node_data() const
    {
        return _middle_node_data;
    }
};

}

#endif /* KATCH_EDGE_INFO_H_ */
