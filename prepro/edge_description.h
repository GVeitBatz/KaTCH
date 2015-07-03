/*
 * katch/datastr/base/edge_description.h
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

#ifndef KATCH_EDGE_DESCRIPTION_H_
#define KATCH_EDGE_DESCRIPTION_H_

namespace katch
{

struct EdgeDescription
{
    using TTF = OrderEdgeData::TTF;

    NodeIterator _src;
    NodeIterator _tgt;
    NodeIterator _middle_node;
    unsigned int _complexity;
    TTF _ttf;
    unsigned int _n_original_edges;

    EdgeDescription(const EdgeDescription&) = delete;
    EdgeDescription& operator= (const EdgeDescription&) = delete;

    EdgeDescription(EdgeDescription&&) = default;
    EdgeDescription& operator= (EdgeDescription&&) = default;

    EdgeDescription()
    : _src(INVALID_NODE_ITERATOR), _tgt(INVALID_NODE_ITERATOR), _middle_node(INVALID_NODE_ITERATOR),
      _complexity(0), _ttf(), _n_original_edges(0)
    {}

    EdgeDescription(const NodeIterator& src, const NodeIterator& tgt,
            const NodeIterator& middle_node, TTF&& ttf, const unsigned int& n_original_edges)
    : _src(src), _tgt(tgt), _middle_node(middle_node),
      _complexity(ttf.size()), _ttf(std::move(ttf)),
      _n_original_edges(n_original_edges)
    {}

    EdgeDescription(const NodeIterator& src, const NodeIterator& tgt,
            const NodeIterator& middle_node, const unsigned int& complexity, const unsigned int& n_original_edges)
    : _src(src), _tgt(tgt), _middle_node(middle_node),
      _complexity(complexity), _ttf(),
      _n_original_edges(n_original_edges)
    {}
};

}

#endif /* KATCH_EDGE_DESCRIPTION_H_ */
