/*
 * Copyright © 2022-2023 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <ios>
#include <cstdlib>
#include <algorithm>

#include "Bytecode.hh"

BytecodeWriter::BytecodeWriter(const filesystem::path &filename)
{
  open(filename, ios::out | ios::binary);
  if (!is_open())
    {
      cerr << R"(Error : Can't open file ")" << filename.string() << R"(" for writing)" << endl;
      exit(EXIT_FAILURE);
    }
}

template<>
BytecodeWriter &
operator<<(BytecodeWriter &code_file, const FCALL_ &instr)
{
  code_file.instructions_positions.push_back(code_file.tellp());

  auto write_member = [&code_file](const auto &member)
  {
    code_file.write(reinterpret_cast<const char *>(&member), sizeof member);
  };

  write_member(instr.op_code);
  write_member(instr.nb_output_arguments);
  write_member(instr.nb_input_arguments);
  write_member(instr.indx);
  write_member(instr.add_input_arguments);
  write_member(instr.row);
  write_member(instr.col);
  write_member(instr.call_type);

  int size = static_cast<int>(instr.func_name.size());
  write_member(size);
  code_file.write(instr.func_name.c_str(), size+1);

  size = static_cast<int>(instr.arg_func_name.size());
  write_member(size);
  code_file.write(instr.arg_func_name.c_str(), size+1);

  return code_file;
}

template<>
BytecodeWriter &
operator<<(BytecodeWriter &code_file, const FBEGINBLOCK_ &instr)
{
  code_file.instructions_positions.push_back(code_file.tellp());

  auto write_member = [&code_file](const auto &member)
  {
    code_file.write(reinterpret_cast<const char *>(&member), sizeof member);
  };

  write_member(instr.op_code);
  write_member(instr.size);
  write_member(instr.type);
  for (int i = 0; i < instr.size; i++)
    {
      write_member(instr.variable[i]);
      write_member(instr.equation[i]);
    }
  if (instr.type == BlockSimulationType::solveTwoBoundariesSimple
      || instr.type == BlockSimulationType::solveTwoBoundariesComplete
      || instr.type == BlockSimulationType::solveBackwardComplete
      || instr.type == BlockSimulationType::solveForwardComplete)
    {
      write_member(instr.is_linear);
      write_member(instr.u_count_int);
    }
  write_member(instr.nb_col_jacob);
  write_member(instr.det_exo_size);
  write_member(instr.exo_size);

  for_each_n(instr.det_exogenous.begin(), instr.det_exo_size, write_member);
  for_each_n(instr.exogenous.begin(), instr.exo_size, write_member);

  return code_file;
}
