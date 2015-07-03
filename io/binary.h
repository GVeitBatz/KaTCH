/*
 * katch/io/binary.h
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

#ifndef KATCH_BINARY_H_
#define KATCH_BINARY_H_

#include <fstream>

namespace katch
{

template <typename T>
void read_value(T& value, std::ifstream& is)
{
    is.read((char*) &value, sizeof(value));
}

uint32_t read_uint32(std::ifstream& is)
{
    uint32_t value;
    read_value(value, is);
    return value;
}

double read_double(std::ifstream& is)
{
    double value;
    read_value(value, is);
    return value;
}

template<typename T>
void write_value(const T value, std::ofstream& os)
{
    os.write((const char*) &value, sizeof(value));
}

void write_uint32(const uint32_t value, std::ofstream& os)
{
    write_value(value, os);
}

void write_double(const double value, std::ofstream& os)
{
    write_value(value, os);
}

}

#endif /* KATCH_BINARY_H_ */
