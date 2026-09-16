/*
 * Copyright © 2026 Dynare Team
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

#include "Utils.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>

#include "Exceptions.hh"

vector<string>
strsplit(string_view str, char delim)
{
  vector<string> result;
  while (true)
    {
      size_t idx {str.find(delim)};
      if (auto sub {str.substr(0, idx)}; !sub.empty())
        result.emplace_back(sub);
      if (idx == string_view::npos)
        break;
      str.remove_prefix(idx + 1);
    }
  return result;
}

filesystem::path
packageDir(const string_view& package)
{
  filesystem::path d;
  for (const auto& it : strsplit(package, '.'))
    d /= "+" + it;
  return d;
}

void
writeToFileIfModified(stringstream& new_contents, const filesystem::path& filename)
{
  ifstream old_file {filename, ios::in | ios::binary};
  if (old_file.is_open()
      && ranges::equal(istreambuf_iterator<char> {old_file}, istreambuf_iterator<char> {},
                       istreambuf_iterator<char> {new_contents}, istreambuf_iterator<char> {}))
    return;
  old_file.close();

  new_contents.seekg(0);

  ofstream new_file {filename, ios::out | ios::binary};
  if (!new_file.is_open())
    throw FileIOException(filename, "writing");
  ranges::copy(istreambuf_iterator<char> {new_contents}, istreambuf_iterator<char> {},
               ostreambuf_iterator<char> {new_file});
  new_file.close();
}

/* Generate a random temporary path, in the current directory. Equivalent to
   boost::filesystem::unique_path(). Both are insecure, but currently there
   is no better portable solution. Maybe in a later C++ standard? */
static filesystem::path
unique_path()
{
  filesystem::path path;
  string possible_characters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  random_device rd;
  mt19937 generator(rd());
  uniform_int_distribution distribution {0, static_cast<int>(possible_characters.size()) - 1};
  do
    {
      constexpr int rand_length = 10;
      string rand_str(rand_length, '\0');
      for (auto& dis : rand_str)
        dis = possible_characters[distribution(generator)];
      path = rand_str;
    }
  while (exists(path));

  return path;
}

void
remove_directory_with_matlab_lock(const filesystem::path& dir)
{
  auto dirStatus {status(dir)};
  if (!exists(dirStatus))
    return;

  if (is_directory(dirStatus))
    for (const auto& e : filesystem::directory_iterator {dir})
      if (e.is_directory())
        remove_directory_with_matlab_lock(e);

  auto tmp {unique_path()};
  rename(dir, tmp);
  remove_all(tmp);
}

void
json_output_with_escape(char c, ostream& output)
{
  switch (c)
    {
    case '\b':
      output << R"(\b)";
      break;

    case '\f':
      output << R"(\f)";
      break;

    case '\n':
      output << R"(\n)";
      break;

    case '\r':
      output << R"(\r)";
      break;

    case '\t':
      output << R"(\t)";
      break;

    case '"':
      output << R"(\")";
      break;

    case '\\':
      output << R"(\\)";
      break;

    default:
      output << c;
      break;
    }
}

void
print_matlab_period(ostream& output, const variant<int, string>& v)
{
  visit([&](const auto& p) { output << p; }, v);
}

void
print_json_period(ostream& output, const variant<int, string>& v)
{
  visit(
      [&]<class T>(const T& p) {
        if constexpr (is_same_v<T, int>)
          output << p;
        else if constexpr (is_same_v<T, string>)
          output << '"' << p << '"';
        else
          static_assert(always_false_v<T>, "Non-exhaustive visitor!");
      },
      v);
}
