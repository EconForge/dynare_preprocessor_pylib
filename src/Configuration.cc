/*
 * Copyright © 2010-2026 Dynare Team
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

#include <array>
#include <fstream>
#include <iostream>
#include <utility>

#ifdef _WIN32
# include <shlobj.h>
#endif

#include "Exceptions.hh"
#include "Configuration.hh"
#include "Utils.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>
#pragma GCC diagnostic pop

Configuration::Path::Path(vector<string> includepath_arg)
{
  if (includepath_arg.empty())
    {
      {
        ostringstream stream;
        stream << "ERROR: The Path must have an Include argument.";
        throw DynareException(stream.str());
      }
    }
  paths["include"] = move(includepath_arg);
}

Configuration::FollowerNode::FollowerNode(string computerName_arg, string port_arg,
                                          int minCpuNbr_arg, int maxCpuNbr_arg, string userName_arg,
                                          string password_arg, string remoteDrive_arg,
                                          string remoteDirectory_arg, string programPath_arg,
                                          string programConfig_arg, string matlabOctavePath_arg,
                                          bool singleCompThread_arg, int numberOfThreadsPerJob_arg,
                                          string operatingSystem_arg) :
    computerName {move(computerName_arg)},
    port {move(port_arg)},
    minCpuNbr {minCpuNbr_arg},
    maxCpuNbr {maxCpuNbr_arg},
    userName {move(userName_arg)},
    password {move(password_arg)},
    remoteDrive {move(remoteDrive_arg)},
    remoteDirectory {move(remoteDirectory_arg)},
    programPath {move(programPath_arg)},
    programConfig {move(programConfig_arg)},
    matlabOctavePath {move(matlabOctavePath_arg)},
    singleCompThread {singleCompThread_arg},
    numberOfThreadsPerJob {numberOfThreadsPerJob_arg},
    operatingSystem {move(operatingSystem_arg)}
{
  if (computerName.empty())
    {
      {
        ostringstream stream;
        stream << "ERROR: The node must have a ComputerName.";
        throw DynareException(stream.str());
      }
    }

  if (!operatingSystem.empty())
    if (operatingSystem != "windows" && operatingSystem != "unix")
      {
        {
        ostringstream stream;
        stream << "ERROR: The OperatingSystem must be either 'unix' or 'windows' (Case Sensitive).";
        throw DynareException(stream.str());
      }
      }
}

Configuration::Cluster::Cluster(member_nodes_t member_nodes_arg) :
    member_nodes {move(member_nodes_arg)}
{
  if (member_nodes.empty())
    {
      {
        ostringstream stream;
        stream << "ERROR: The cluster must have at least one member node.";
        throw DynareException(stream.str());
      }
    }
}

Configuration::Configuration(bool parallel_arg, bool parallel_test_arg,
                             bool parallel_follower_open_mode_arg, bool parallel_use_psexec_arg,
                             string cluster_name_arg) :
    parallel {parallel_arg},
    parallel_test {parallel_test_arg},
    parallel_follower_open_mode {parallel_follower_open_mode_arg},
    parallel_use_psexec {parallel_use_psexec_arg},
    cluster_name {move(cluster_name_arg)}
{
}

void
Configuration::getConfigFileInfo(const filesystem::path& conffile_option,
                                 WarningConsolidation& warnings)
{
  using namespace boost;

  filesystem::path config_file {conffile_option};

  if (config_file.empty())
    {
      config_file = findConfigFile("dynare.ini");

      if (config_file.empty()) // Try old default location (Dynare ⩽ 5) for backward compatibility
        {
          filesystem::path old_default_config_file;
#ifdef _WIN32
          array<wchar_t, MAX_PATH + 1> appdata;
          if (SHGetFolderPathW(nullptr, CSIDL_APPDATA | CSIDL_FLAG_DONT_VERIFY, nullptr,
                               SHGFP_TYPE_CURRENT, appdata.data())
              == S_OK)
            old_default_config_file = filesystem::path {appdata.data()} / "dynare.ini";
#else
          if (auto home = getenv("HOME"); home)
            old_default_config_file = filesystem::path {home} / ".dynare";
#endif
          if (!old_default_config_file.empty() && exists(old_default_config_file))
            {
              warnings << "WARNING: the location " << old_default_config_file.string()
                       << " for the configuration file is obsolete; please see the reference"
                       << " manual for the new location." << '\n';
              config_file = old_default_config_file;
            }
        }
    }

  if (config_file.empty())
    {
      if (parallel || parallel_test)
        {
          {
        ostringstream stream;
        stream << "ERROR: the parallel or parallel_test option was passed but no configuration "
               << "file was found";
        throw DynareException(stream.str());
      }
        }
      else
        return;
    }

  ifstream configFile;
  configFile.open(config_file, fstream::in);

  if (!configFile.is_open())
    {
      {
        ostringstream stream;
        stream << "ERROR: Couldn't open configuration file " << config_file.string();
        throw DynareException(stream.str());
      }
    }

  string name, computerName, port, userName, password, remoteDrive, remoteDirectory, programPath,
      programConfig, matlabOctavePath, operatingSystem;
  vector<string> includepath;
  int minCpuNbr {0}, maxCpuNbr {0};
  int numberOfThreadsPerJob {1};
  bool singleCompThread {false};
  member_nodes_t member_nodes;

  bool inHooks {false}, inNode {false}, inCluster {false}, inPaths {false};

  while (configFile.good())
    {
      string line;
      getline(configFile, line);
      trim(line);
      if (line.empty() || !line.compare(0, 1, "#"))
        continue;

      if (line == "[node]" || line == "[cluster]" || line == "[hooks]" || line == "[paths]")
        {
          if (!includepath.empty())
            // we were just in [paths]
            addPathsConfFileElement(includepath);
          else
            // we were just in [node] or [cluster]
            addParallelConfFileElement(
                inNode, inCluster, member_nodes, name, computerName, port, minCpuNbr, maxCpuNbr,
                userName, password, remoteDrive, remoteDirectory, programPath, programConfig,
                matlabOctavePath, singleCompThread, numberOfThreadsPerJob, operatingSystem);

          //! Reset communication vars / option defaults
          if (line == "[hooks]")
            {
              inHooks = true;
              inNode = false;
              inCluster = false;
              inPaths = false;
            }
          else if (line == "[node]")
            {
              inHooks = false;
              inNode = true;
              inCluster = false;
              inPaths = false;
            }
          else if (line == "[paths]")
            {
              inHooks = false;
              inNode = false;
              inCluster = false;
              inPaths = true;
            }
          else
            {
              inHooks = false;
              inNode = false;
              inCluster = true;
              inPaths = false;
            }

          name = userName = computerName = port = password = remoteDrive = remoteDirectory
              = programPath = programConfig = matlabOctavePath = operatingSystem = "";
          includepath.clear();
          minCpuNbr = maxCpuNbr = 0;
          numberOfThreadsPerJob = 1;
          singleCompThread = false;
          member_nodes.clear();
        }
      else
        {
          vector<string> tokenizedLine;
          split(tokenizedLine, line, is_any_of("="));
          if (tokenizedLine.size() != 2)
            {
              {
        ostringstream stream;
        stream << "ERROR (in config file): Options should be formatted as 'option = value'.";
        throw DynareException(stream.str());
      }
            }
          trim(tokenizedLine.front());
          trim(tokenizedLine.back());

          if (inHooks)
            if (tokenizedLine.front() == "GlobalInitFile")
              if (global_init_file.empty())
                global_init_file = tokenizedLine.back();
              else
                {
                  {
        ostringstream stream;
        stream << "ERROR: May not have more than one GlobalInitFile option in [hooks] block.";
        throw DynareException(stream.str());
      }
                }
            else
              {
                {
        ostringstream stream;
        stream << "ERROR: Unrecognized option " << tokenizedLine.front()
                     << " in [hooks] block.";
        throw DynareException(stream.str());
      }
              }
          else if (inPaths)
            if (tokenizedLine.front() == "Include")
              if (includepath.empty())
                {
                  vector<string> tokenizedPath;
                  split(tokenizedPath, tokenizedLine.back(), is_any_of(":"), token_compress_on);
                  for (auto& it : tokenizedPath)
                    if (!it.empty())
                      {
                        trim(it);
                        includepath.push_back(it);
                      }
                }
              else
                {
                  {
        ostringstream stream;
        stream << "ERROR: May not have more than one Include option in [paths] block.";
        throw DynareException(stream.str());
      }
                }
            else
              {
                {
        ostringstream stream;
        stream << "ERROR: Unrecognized option " << tokenizedLine.front()
                     << " in [paths] block.";
        throw DynareException(stream.str());
      }
              }
          else if (tokenizedLine.front() == "Name")
            name = tokenizedLine.back();
          else if (tokenizedLine.front() == "CPUnbr")
            {
              vector<string> tokenizedCpuNbr;
              split(tokenizedCpuNbr, tokenizedLine.back(), is_any_of(":"));
              try
                {
                  if (tokenizedCpuNbr.size() == 1)
                    {
                      minCpuNbr = 1;
                      maxCpuNbr = stoi(tokenizedCpuNbr.front());
                    }
                  else if (tokenizedCpuNbr.size() == 2 && tokenizedCpuNbr[0].at(0) == '['
                           && tokenizedCpuNbr[1].at(tokenizedCpuNbr[1].size() - 1) == ']')
                    {
                      tokenizedCpuNbr[0].erase(0, 1);
                      tokenizedCpuNbr[1].erase(tokenizedCpuNbr[1].size() - 1, 1);
                      minCpuNbr = stoi(tokenizedCpuNbr[0]);
                      maxCpuNbr = stoi(tokenizedCpuNbr[1]);
                    }
                }
              catch (const invalid_argument&)
                {
                  {
        ostringstream stream;
        stream << "ERROR: Could not convert value to integer for CPUnbr.";
        throw DynareException(stream.str());
      }
                }

              if (minCpuNbr <= 0 || maxCpuNbr <= 0)
                {
                  {
        ostringstream stream;
        stream << "ERROR: Syntax for the CPUnbr option is as follows:" << '\n'
                       << "       1) CPUnbr = <int>" << '\n'
                       << "    or 2) CPUnbr = [<int>:<int>]" << '\n'
                       << "       where <int> is an Integer > 0.";
        throw DynareException(stream.str());
      }
                }

              minCpuNbr--;
              maxCpuNbr--;
              if (minCpuNbr > maxCpuNbr)
                {
                  int tmp = maxCpuNbr;
                  maxCpuNbr = minCpuNbr;
                  minCpuNbr = tmp;
                }
            }
          else if (tokenizedLine.front() == "Port")
            port = tokenizedLine.back();
          else if (tokenizedLine.front() == "ComputerName")
            computerName = tokenizedLine.back();
          else if (tokenizedLine.front() == "UserName")
            userName = tokenizedLine.back();
          else if (tokenizedLine.front() == "Password")
            password = tokenizedLine.back();
          else if (tokenizedLine.front() == "RemoteDrive")
            remoteDrive = tokenizedLine.back();
          else if (tokenizedLine.front() == "RemoteDirectory")
            remoteDirectory = tokenizedLine.back();
          else if (tokenizedLine.front() == "DynarePath" || tokenizedLine.front() == "ProgramPath")
            programPath = tokenizedLine.back();
          else if (tokenizedLine.front() == "ProgramConfig")
            programConfig = tokenizedLine.back();
          else if (tokenizedLine.front() == "MatlabOctavePath")
            matlabOctavePath = tokenizedLine.back();
          else if (tokenizedLine.front() == "NumberOfThreadsPerJob")
            numberOfThreadsPerJob = stoi(tokenizedLine.back());
          else if (tokenizedLine.front() == "SingleCompThread")
            if (tokenizedLine.back() == "true")
              singleCompThread = true;
            else if (tokenizedLine.back() == "false")
              singleCompThread = false;
            else
              {
                {
        ostringstream stream;
        stream << "ERROR (in config file): The value passed to SingleCompThread may only be "
                        "'true' or 'false'.";
        throw DynareException(stream.str());
      }
              }
          else if (tokenizedLine.front() == "OperatingSystem")
            operatingSystem = tokenizedLine.back();
          else if (tokenizedLine.front() == "Members")
            {
              char_separator sep(" ,;", "()", drop_empty_tokens);
              tokenizer tokens(tokenizedLine.back(), sep);
              string node_name;
              for (bool begin_weight {false}; const auto& token : tokens)
                {
                  if (token == "(")
                    {
                      begin_weight = true;
                      continue;
                    }
                  else if (token == ")")
                    {
                      node_name.clear();
                      begin_weight = false;
                      continue;
                    }

                  if (!begin_weight)
                    {
                      if (!node_name.empty())
                        {
                          if (member_nodes.contains(node_name))
                            {
                              {
        ostringstream stream;
        stream << "ERROR (in config file): Node entered twice in specification "
                                      "of cluster.";
        throw DynareException(stream.str());
      }
                            }
                          else
                            member_nodes[node_name] = 1.0;
                        }
                      node_name = token;
                    }
                  else
                    try
                      {
                        auto weight = stod(token);
                        if (weight <= 0)
                          {
                            {
        ostringstream stream;
        stream << "ERROR (in config file): Misspecification of weights passed to "
                                    "Members option.";
        throw DynareException(stream.str());
      }
                          }
                        member_nodes[node_name] = weight;
                      }
                    catch (const invalid_argument&)
                      {
                        {
        ostringstream stream;
        stream << "ERROR (in config file): Misspecification of weights passed to "
                                "Members option.";
        throw DynareException(stream.str());
      }
                      }
                }
              if (!node_name.empty())
                {
                  if (!member_nodes.contains(node_name))
                    member_nodes[node_name] = 1.0;
                  else
                    {
                      {
        ostringstream stream;
        stream << "ERROR (in config file): Node entered twice in specification of "
                              "cluster.";
        throw DynareException(stream.str());
      }
                    }
                }
            }
          else
            {
              {
        ostringstream stream;
        stream << "ERROR (in config file): Option " << tokenizedLine.front() << " is invalid.";
        throw DynareException(stream.str());
      }
            }
        }
    }

  if (!includepath.empty())
    addPathsConfFileElement(includepath);
  else
    addParallelConfFileElement(inNode, inCluster, member_nodes, name, computerName, port, minCpuNbr,
                               maxCpuNbr, userName, password, remoteDrive, remoteDirectory,
                               programPath, programConfig, matlabOctavePath, singleCompThread,
                               numberOfThreadsPerJob, operatingSystem);

  configFile.close();
}

void
Configuration::addPathsConfFileElement(vector<string> includepath)
{
  if (includepath.empty())
    {
      {
        ostringstream stream;
        stream << "ERROR: The path to be included must be passed to the Include option.";
        throw DynareException(stream.str());
      }
    }
  else
    paths.emplace_back(move(includepath));
}

void
Configuration::addParallelConfFileElement(bool inNode, bool inCluster,
                                          const member_nodes_t& member_nodes, const string& name,
                                          const string& computerName, const string& port,
                                          int minCpuNbr, int maxCpuNbr, const string& userName,
                                          const string& password, const string& remoteDrive,
                                          const string& remoteDirectory, const string& programPath,
                                          const string& programConfig,
                                          const string& matlabOctavePath, bool singleCompThread,
                                          int numberOfThreadsPerJob, const string& operatingSystem)
{
  //! ADD NODE
  if (inNode)
    if (!member_nodes.empty())
      {
        {
        ostringstream stream;
        stream << "Invalid option passed to [node].";
        throw DynareException(stream.str());
      }
      }
    else if (name.empty() || follower_nodes.contains(name))
      {
        {
        ostringstream stream;
        stream << "ERROR: Every node must be assigned a unique name.";
        throw DynareException(stream.str());
      }
      }
    else
      follower_nodes.try_emplace(name, computerName, port, minCpuNbr, maxCpuNbr, userName, password,
                                 remoteDrive, remoteDirectory, programPath, programConfig,
                                 matlabOctavePath, singleCompThread, numberOfThreadsPerJob,
                                 operatingSystem);
  //! ADD CLUSTER
  else if (inCluster)
    {
      if (minCpuNbr > 0 || maxCpuNbr > 0 || !userName.empty() || !password.empty()
          || !remoteDrive.empty() || !remoteDirectory.empty() || !programPath.empty()
          || !programConfig.empty() || !matlabOctavePath.empty() || !operatingSystem.empty())
        {
          {
        ostringstream stream;
        stream << "Invalid option passed to [cluster].";
        throw DynareException(stream.str());
      }
        }
      else if (name.empty() || clusters.contains(name))
        {
          {
        ostringstream stream;
        stream << "ERROR: The cluster must be assigned a unique name.";
        throw DynareException(stream.str());
      }
        }
      else
        {
          if (clusters.empty())
            firstClusterName = name;
          clusters.emplace(name, member_nodes);
        }
    }
}

void
Configuration::checkPass([[maybe_unused]] WarningConsolidation& warnings) const
{
  if (!parallel && !parallel_test)
    return;

  //! Check Follower Nodes
  if (follower_nodes.empty())
    {
      {
        ostringstream stream;
        stream << "ERROR: At least one node must be defined in the config file.";
        throw DynareException(stream.str());
      }
    }

  for (const auto& follower_node : follower_nodes)
    {
#if !defined(_WIN32) && !defined(__CYGWIN32__)
      // For Linux/Mac, check that cpuNbr starts at 0
      if (follower_node.second.minCpuNbr != 0)
        warnings << "WARNING: On Unix-based operating systems, you cannot specify the CPU that is "
                 << "used in parallel processing. This will be adjusted for you such that the "
                 << "same number of CPUs are used." << '\n';
#endif
      if (!follower_node.second.port.empty())
        try
          {
            stoi(follower_node.second.port);
          }
        catch (const invalid_argument&)
          {
            {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first << "): the port must be an integer.";
        throw DynareException(stream.str());
      }
          }
      if (follower_node.second.computerName == "localhost") // We are working locally
        {
          if (!follower_node.second.remoteDrive.empty())
            {
              {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first
                   << "): the RemoteDrive option may not be passed for a local node.";
        throw DynareException(stream.str());
      }
            }
          if (!follower_node.second.remoteDirectory.empty())
            {
              {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first
                   << "): the RemoteDirectory option may not be passed for a local node.";
        throw DynareException(stream.str());
      }
            }
        }
      else
        {
          if (follower_node.second.userName.empty())
            {
              {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first
                   << "): the UserName option must be passed for every remote node.";
        throw DynareException(stream.str());
      }
            }
          if (follower_node.second.operatingSystem == "windows")
            {
              if (follower_node.second.password.empty())
                {
                  {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first
                       << "): the Password option must be passed under Windows for every remote "
                          "node.";
        throw DynareException(stream.str());
      }
                }
              if (follower_node.second.remoteDrive.empty())
                {
                  {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first
                       << "): the RemoteDrive option must be passed under Windows for every remote "
                          "node.";
        throw DynareException(stream.str());
      }
                }
            }
#if defined(_WIN32) || defined(__CYGWIN32__)
          if (follower_node.second.operatingSystem.empty())
            {
              if (follower_node.second.password.empty())
                {
                  {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first
                       << "): the Password option must be passed under Windows for every remote "
                          "node.";
        throw DynareException(stream.str());
      }
                }
              if (follower_node.second.remoteDrive.empty())
                {
                  {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first
                       << "): the RemoteDrive option must be passed under Windows for every remote "
                          "node.";
        throw DynareException(stream.str());
      }
                }
            }
#endif
          if (follower_node.second.remoteDirectory.empty())
            {
              {
        ostringstream stream;
        stream << "ERROR (node " << follower_node.first
                   << "): the RemoteDirectory must be specified for every remote node.";
        throw DynareException(stream.str());
      }
            }
        }
    }

  //! Check Clusters
  if (clusters.empty())
    {
      {
        ostringstream stream;
        stream << "ERROR: At least one cluster must be defined in the config file.";
        throw DynareException(stream.str());
      }
    }

  if (!cluster_name.empty() && !clusters.contains(cluster_name))
    {
      {
        ostringstream stream;
        stream << "ERROR: Cluster Name " << cluster_name << " was not found in the config file.";
        throw DynareException(stream.str());
      }
    }

  for (const auto& cluster : clusters)
    for (const auto& itmn : cluster.second.member_nodes)
      if (!follower_nodes.contains(itmn.first))
        {
          {
        ostringstream stream;
        stream << "Error: node " << itmn.first << " specified in cluster " << cluster.first
               << " was not found";
        throw DynareException(stream.str());
      }
        }
}

void
Configuration::transformPass()
{
  /* If the user did not specify the GlobalInitFile option, use global_init.m in configuration
     directory if it exists */
  if (auto default_global_init_file = findConfigFile("global_init.m");
      global_init_file.empty() && !default_global_init_file.empty())
    global_init_file = default_global_init_file.string();

  if (!parallel && !parallel_test)
    return;

#if !defined(_WIN32) && !defined(__CYGWIN32__)
  // For Linux/Mac, check that cpuNbr starts at 0
  for (auto& it : follower_nodes)
    if (it.second.minCpuNbr != 0)
      {
        it.second.maxCpuNbr = it.second.maxCpuNbr - it.second.minCpuNbr;
        it.second.minCpuNbr = 0;
      }
#endif

  auto& cluster = cluster_name.empty() ? clusters.at(firstClusterName) : clusters.at(cluster_name);

  double weight_denominator {0.0};
  for (const auto& [name, weight] : cluster.member_nodes)
    weight_denominator += weight;

  for (auto& [name, weight] : cluster.member_nodes)
    weight /= weight_denominator;
}

vector<filesystem::path>
Configuration::getIncludePaths() const
{
  vector<filesystem::path> include_paths;
  for (const auto& path : paths)
    for (const auto& mapit : path.get_paths())
      for (const auto& vecit : mapit.second)
        include_paths.emplace_back(vecit);
  return include_paths;
}

void
Configuration::writeHooks(ostream& output) const
{
  if (!global_init_file.empty())
    output << "options_.global_init_file = '" << global_init_file << "';" << '\n';
}

void
Configuration::writeCluster(ostream& output) const
{
  if (!parallel && !parallel_test)
    return;

  const auto& cluster
      = cluster_name.empty() ? clusters.at(firstClusterName) : clusters.at(cluster_name);

  for (int i {1}; const auto& [name, node] : follower_nodes)
    {
      if (!cluster.member_nodes.contains(name))
        continue; // Skip nodes not in the selected cluster

      output << "options_.parallel";
      if (i > 1)
        output << "(" << i << ")";
      i++;
      output << " = struct('Local', " << noboolalpha << (node.computerName == "localhost") << ", "
             << "'ComputerName', '" << node.computerName << "', "
             << "'Port', '" << node.port << "', "
             << "'CPUnbr', [" << node.minCpuNbr << ":" << node.maxCpuNbr << "], "
             << "'UserName', '" << node.userName << "', "
             << "'Password', '" << node.password << "', "
             << "'RemoteDrive', '" << node.remoteDrive << "', "
             << "'RemoteDirectory', '" << node.remoteDirectory
             << "', "
             // The following should be switched back to “ProgramPath” once we move to Dragonfly
             << "'DynarePath', '" << node.programPath << "', "
             << "'ProgramConfig', '" << node.programConfig << "', "
             << "'MatlabOctavePath', '" << node.matlabOctavePath << "', "
             << "'OperatingSystem', '" << node.operatingSystem << "', "
             << "'NodeWeight', '" << cluster.member_nodes.at(name) << "', "
             << "'NumberOfThreadsPerJob', " << node.numberOfThreadsPerJob << ", "
             << "'SingleCompThread', '" << boolalpha << node.singleCompThread << "');" << '\n';
    }

  // Default values for the following two are both in DynareMain.cc and
  // matlab/default_option_values.m
  if (parallel_follower_open_mode)
    output << "options_.parallel_info.leaveSlaveOpen = 1;" << '\n';
  if (!parallel_use_psexec)
    output << "options_.parallel_info.use_psexec = false;" << '\n';

  output << "options_.parallel_info.console_mode= isoctave;" << '\n';

  output << "InitializeComputationalEnvironment();" << '\n';
  if (parallel_test)
    output
        << "ErrorCode = AnalyseComputationalEnvironment(options_.parallel, options_.parallel_info);"
        << '\n'
        << "disp(['AnalyseComputationalEnvironment returned with Error Code: ' "
           "num2str(ErrorCode)]);"
        << '\n'
        << "diary off;" << '\n'
        << "return;" << '\n';
}

void
Configuration::writeEndParallel(ostream& output) const
{
  if ((!parallel && !parallel_test) || !parallel_follower_open_mode)
    return;

  output << "if options_.parallel_info.leaveSlaveOpen == 1" << '\n'
         << "     closeSlave(options_.parallel,options_.parallel_info.RemoteTmpFolder);" << '\n'
         << "end" << '\n';
}

filesystem::path
Configuration::findConfigFile(const string& filename)
{
#ifdef _WIN32
  array<wchar_t, MAX_PATH + 1> appdata;
  if (SHGetFolderPathW(nullptr, CSIDL_APPDATA | CSIDL_FLAG_DONT_VERIFY, nullptr, SHGFP_TYPE_CURRENT,
                       appdata.data())
      == S_OK)
    {
      filesystem::path candidate {filesystem::path {appdata.data()} / "dynare" / filename};
      if (exists(candidate))
        return candidate;
    }
#else
  filesystem::path xdg_config_home;
  if (auto xdg_config_home_env = getenv("XDG_CONFIG_HOME"); xdg_config_home_env)
    xdg_config_home = xdg_config_home_env;
  if (auto home = getenv("HOME"); xdg_config_home.empty() && home)
    xdg_config_home = filesystem::path {home} / ".config";

  if (!xdg_config_home.empty())
    {
      filesystem::path candidate {xdg_config_home / "dynare" / filename};
      if (exists(candidate))
        return candidate;
    }

  string xdg_config_dirs;
  if (auto xdg_config_dirs_env = getenv("XDG_CONFIG_DIRS"); xdg_config_dirs_env)
    xdg_config_dirs = xdg_config_dirs_env;
  if (xdg_config_dirs.empty())
    xdg_config_dirs = "/etc/xdg";
  for (const auto& dir : strsplit(xdg_config_dirs, ':'))
    {
      filesystem::path candidate {filesystem::path {dir} / "dynare" / filename};
      if (exists(candidate))
        return candidate;
    }
#endif

  return {};
}
