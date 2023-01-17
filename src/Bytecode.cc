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

  code_file.write(reinterpret_cast<const char *>(&instr.op_code), sizeof instr.op_code);
  code_file.write(reinterpret_cast<const char *>(&instr.nb_output_arguments), sizeof instr.nb_output_arguments);
  code_file.write(reinterpret_cast<const char *>(&instr.nb_input_arguments), sizeof instr.nb_input_arguments);
  code_file.write(reinterpret_cast<const char *>(&instr.indx), sizeof instr.indx);
  code_file.write(reinterpret_cast<const char *>(&instr.add_input_arguments), sizeof instr.add_input_arguments);
  code_file.write(reinterpret_cast<const char *>(&instr.row), sizeof instr.row);
  code_file.write(reinterpret_cast<const char *>(&instr.col), sizeof instr.col);
  code_file.write(reinterpret_cast<const char *>(&instr.call_type), sizeof instr.call_type);

  int size = static_cast<int>(instr.func_name.size());
  code_file.write(reinterpret_cast<char *>(&size), sizeof size);
  const char *name = instr.func_name.c_str();
  code_file.write(name, size);

  size = static_cast<int>(instr.arg_func_name.size());
  code_file.write(reinterpret_cast<char *>(&size), sizeof size);
  name = instr.arg_func_name.c_str();
  code_file.write(name, size);

  return code_file;
}

template<>
BytecodeWriter &
operator<<(BytecodeWriter &code_file, const FBEGINBLOCK_ &instr)
{
  code_file.instructions_positions.push_back(code_file.tellp());

  code_file.write(reinterpret_cast<const char *>(&instr.op_code), sizeof instr.op_code);
  code_file.write(reinterpret_cast<const char *>(&instr.size), sizeof instr.size);
  code_file.write(reinterpret_cast<const char *>(&instr.type), sizeof instr.type);
  for (int i = 0; i < instr.size; i++)
    {
      code_file.write(reinterpret_cast<const char *>(&instr.variable[i]), sizeof instr.variable[i]);
      code_file.write(reinterpret_cast<const char *>(&instr.equation[i]), sizeof instr.equation[i]);
    }
  if (instr.type == BlockSimulationType::solveTwoBoundariesSimple
      || instr.type == BlockSimulationType::solveTwoBoundariesComplete
      || instr.type == BlockSimulationType::solveBackwardComplete
      || instr.type == BlockSimulationType::solveForwardComplete)
    {
      code_file.write(reinterpret_cast<const char *>(&instr.is_linear), sizeof instr.is_linear);
      code_file.write(reinterpret_cast<const char *>(&instr.endo_nbr), sizeof instr.endo_nbr);
      code_file.write(reinterpret_cast<const char *>(&instr.u_count_int), sizeof instr.u_count_int);
    }
  code_file.write(reinterpret_cast<const char *>(&instr.nb_col_jacob), sizeof instr.nb_col_jacob);
  code_file.write(reinterpret_cast<const char *>(&instr.det_exo_size), sizeof instr.det_exo_size);
  code_file.write(reinterpret_cast<const char *>(&instr.exo_size), sizeof instr.exo_size);

  for (int i{0}; i < instr.det_exo_size; i++)
    code_file.write(reinterpret_cast<const char *>(&instr.det_exogenous[i]), sizeof instr.det_exogenous[i]);
  for (int i{0}; i < instr.exo_size; i++)
    code_file.write(reinterpret_cast<const char *>(&instr.exogenous[i]), sizeof instr.exogenous[i]);

  return code_file;
}
