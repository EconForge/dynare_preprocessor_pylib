/*
 * Copyright © 2010-2023 Dynare Team
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
#include <fstream>
#include <utility>
#include <vector>

#include "ConfigFile.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>
#pragma GCC diagnostic pop

using namespace std;

Hook::Hook(string global_init_file_arg)
{
  if (global_init_file_arg.empty())
    {
      cerr << "ERROR: The Hook must have a Global Initialization File argument." << endl;
      exit(EXIT_FAILURE);
    }
  hooks["global_init_file"] = move(global_init_file_arg);
}

Path::Path(vector<string> includepath_arg)
{
  if (includepath_arg.empty())
    {
      cerr << "ERROR: The Path must have an Include argument." << endl;
      exit(EXIT_FAILURE);
    }
  paths["include"] = move(includepath_arg);
}

FollowerNode::FollowerNode(string computerName_arg, string port_arg, int minCpuNbr_arg, int maxCpuNbr_arg, string userName_arg,
                           string password_arg, string remoteDrive_arg, string remoteDirectory_arg,
                           string programPath_arg, string programConfig_arg, string matlabOctavePath_arg, bool singleCompThread_arg,
                           int numberOfThreadsPerJob_arg, string operatingSystem_arg) :
  computerName{move(computerName_arg)},
  port{move(port_arg)},
  minCpuNbr{minCpuNbr_arg},
  maxCpuNbr{maxCpuNbr_arg},
  userName{move(userName_arg)},
  password{move(password_arg)},
  remoteDrive{move(remoteDrive_arg)},
  remoteDirectory{move(remoteDirectory_arg)},
  programPath{move(programPath_arg)},
  programConfig{move(programConfig_arg)},
  matlabOctavePath{move(matlabOctavePath_arg)},
  singleCompThread{singleCompThread_arg},
  numberOfThreadsPerJob{numberOfThreadsPerJob_arg},
  operatingSystem{move(operatingSystem_arg)}
{
  if (computerName.empty())
    {
      cerr << "ERROR: The node must have a ComputerName." << endl;
      exit(EXIT_FAILURE);
    }

  if (!operatingSystem.empty())
    if (operatingSystem != "windows" && operatingSystem != "unix")
      {
        cerr << "ERROR: The OperatingSystem must be either 'unix' or 'windows' (Case Sensitive)." << endl;
        exit(EXIT_FAILURE);
      }
}

Cluster::Cluster(member_nodes_t member_nodes_arg) :
  member_nodes{move(member_nodes_arg)}
{
  if (member_nodes.empty())
    {
      cerr << "ERROR: The cluster must have at least one member node." << endl;
      exit(EXIT_FAILURE);
    }
}

ConfigFile::ConfigFile(bool parallel_arg, bool parallel_test_arg,
                       bool parallel_follower_open_mode_arg, bool parallel_use_psexec_arg,
                       string cluster_name_arg) :
  parallel{parallel_arg}, parallel_test{parallel_test_arg},
  parallel_follower_open_mode{parallel_follower_open_mode_arg},
  parallel_use_psexec{parallel_use_psexec_arg},
  cluster_name{move(cluster_name_arg)}
{
}

void
ConfigFile::getConfigFileInfo(const filesystem::path &config_file)
{
  using namespace boost;
  ifstream configFile;

  if (config_file.empty())
    {
      filesystem::path defaultConfigFile;
      // Test OS and try to open default file
#if defined(_WIN32) || defined(__CYGWIN32__)
      if (auto appdata = getenv("APPDATA");
          appdata)
        defaultConfigFile = filesystem::path{appdata} / "dynare.ini";
      else
        {
          if (parallel || parallel_test)
            cerr << "ERROR: ";
          else
            cerr << "WARNING: ";
          cerr << "APPDATA environment variable not found." << endl;

          if (parallel || parallel_test)
            exit(EXIT_FAILURE);
        }
#else
      if (auto home = getenv("HOME");
          home)
        defaultConfigFile = filesystem::path{home} / ".dynare";
      else
        {
          if (parallel || parallel_test)
            cerr << "ERROR: ";
          else
            cerr << "WARNING: ";
          cerr << "HOME environment variable not found." << endl;
          if (parallel || parallel_test)
            exit(EXIT_FAILURE);
        }
#endif
      configFile.open(defaultConfigFile, fstream::in);
      if (!configFile.is_open())
        if (parallel || parallel_test)
          {
            cerr << "ERROR: Could not open the default config file (" << defaultConfigFile.string() << ")" << endl;
            exit(EXIT_FAILURE);
          }
        else
          return;
    }
  else
    {
      configFile.open(config_file, fstream::in);
      if (!configFile.is_open())
        {
          cerr << "ERROR: Couldn't open file " << config_file.string() << endl;;
          exit(EXIT_FAILURE);
        }
    }

  string name, computerName, port, userName, password, remoteDrive,
    remoteDirectory, programPath, programConfig, matlabOctavePath,
    operatingSystem, global_init_file;
  vector<string> includepath;
  int minCpuNbr{0}, maxCpuNbr{0};
  int numberOfThreadsPerJob{1};
  bool singleCompThread{false};
  member_nodes_t member_nodes;

  bool inHooks{false}, inNode{false}, inCluster{false}, inPaths{false};

  while (configFile.good())
    {
      string line;
      getline(configFile, line);
      trim(line);
      if (line.empty() || !line.compare(0, 1, "#"))
        continue;

      if (line == "[node]"
          || line == "[cluster]"
          || line == "[hooks]"
          || line == "[paths]")
        {
          if (!global_init_file.empty())
            // we were just in [hooks]
            addHooksConfFileElement(global_init_file);
          else if (!includepath.empty())
            // we were just in [paths]
            addPathsConfFileElement(includepath);
          else
            // we were just in [node] or [cluster]
            addParallelConfFileElement(inNode, inCluster, member_nodes, name,
                                       computerName, port, minCpuNbr, maxCpuNbr, userName,
                                       password, remoteDrive, remoteDirectory,
                                       programPath, programConfig, matlabOctavePath, singleCompThread,
                                       numberOfThreadsPerJob, operatingSystem);

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

          name = userName = computerName = port = password = remoteDrive
            = remoteDirectory = programPath = programConfig = matlabOctavePath
            = operatingSystem = global_init_file = "";
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
              cerr << "ERROR (in config file): Options should be formatted as 'option = value'." << endl;
              exit(EXIT_FAILURE);
            }
          trim(tokenizedLine.front());
          trim(tokenizedLine.back());

          if (inHooks)
            if (tokenizedLine.front() == "GlobalInitFile")
              if (global_init_file.empty())
                global_init_file = tokenizedLine.back();
              else
                {
                  cerr << "ERROR: May not have more than one GlobalInitFile option in [hooks] block." << endl;
                  exit(EXIT_FAILURE);
                }
            else
              {
                cerr << "ERROR: Unrecognized option " << tokenizedLine.front() << " in [hooks] block." << endl;
                exit(EXIT_FAILURE);
              }
          else if (inPaths)
            if (tokenizedLine.front() == "Include")
              if (includepath.empty())
                {
                  vector<string> tokenizedPath;
                  split(tokenizedPath, tokenizedLine.back(), is_any_of(":"), token_compress_on);
                  for (auto &it : tokenizedPath)
                    if (!it.empty())
                      {
                        trim(it);
                        includepath.push_back(it);
                      }
                }
              else
                {
                  cerr << "ERROR: May not have more than one Include option in [paths] block." << endl;
                  exit(EXIT_FAILURE);
                }
            else
              {
                cerr << "ERROR: Unrecognized option " << tokenizedLine.front() << " in [paths] block." << endl;
                exit(EXIT_FAILURE);
              }
          else
            if (tokenizedLine.front() == "Name")
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
                    else if (tokenizedCpuNbr.size() == 2
                             && tokenizedCpuNbr[0].at(0) == '['
                             && tokenizedCpuNbr[1].at(tokenizedCpuNbr[1].size()-1) == ']')
                      {
                        tokenizedCpuNbr[0].erase(0, 1);
                        tokenizedCpuNbr[1].erase(tokenizedCpuNbr[1].size()-1, 1);
                        minCpuNbr = stoi(tokenizedCpuNbr[0]);
                        maxCpuNbr = stoi(tokenizedCpuNbr[1]);
                      }
                  }
                catch (const invalid_argument &)
                  {
                    cerr << "ERROR: Could not convert value to integer for CPUnbr." << endl;
                    exit(EXIT_FAILURE);
                  }

                if (minCpuNbr <= 0 || maxCpuNbr <= 0)
                  {
                    cerr << "ERROR: Syntax for the CPUnbr option is as follows:" << endl
                         << "       1) CPUnbr = <int>" << endl
                         << "    or 2) CPUnbr = [<int>:<int>]" << endl
                         << "       where <int> is an Integer > 0." << endl;
                    exit(EXIT_FAILURE);
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
            else if (tokenizedLine.front() == "DynarePath"
                     || tokenizedLine.front() == "ProgramPath")
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
                  cerr << "ERROR (in config file): The value passed to SingleCompThread may only be 'true' or 'false'." << endl;
                  exit(EXIT_FAILURE);
                }
            else if (tokenizedLine.front() == "OperatingSystem")
              operatingSystem = tokenizedLine.back();
            else if (tokenizedLine.front() == "Members")
              {
                char_separator sep(" ,;", "()", drop_empty_tokens);
                tokenizer tokens(tokenizedLine.back(), sep);
                string node_name;
                for (bool begin_weight{false};
                     const auto &token : tokens)
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
                          if (member_nodes.contains(node_name))
                            {
                              cerr << "ERROR (in config file): Node entered twice in specification of cluster." << endl;
                              exit(EXIT_FAILURE);
                            }
                          else
                            member_nodes[node_name] = 1.0;
                        node_name = token;
                      }
                    else
                      try
                        {
                          auto weight = stod(token);
                          if (weight <= 0)
                            {
                              cerr << "ERROR (in config file): Misspecification of weights passed to Members option." << endl;
                              exit(EXIT_FAILURE);
                            }
                          member_nodes[node_name] = weight;
                        }
                      catch (const invalid_argument &)
                        {
                          cerr << "ERROR (in config file): Misspecification of weights passed to Members option." << endl;
                          exit(EXIT_FAILURE);
                        }
                  }
                if (!node_name.empty())
                  if (!member_nodes.contains(node_name))
                    member_nodes[node_name] = 1.0;
                  else
                    {
                      cerr << "ERROR (in config file): Node entered twice in specification of cluster." << endl;
                      exit(EXIT_FAILURE);
                    }
              }
            else
              {
                cerr << "ERROR (in config file): Option " << tokenizedLine.front() << " is invalid." << endl;
                exit(EXIT_FAILURE);
              }
        }
    }

  if (!global_init_file.empty())
    addHooksConfFileElement(global_init_file);
  else if (!includepath.empty())
    addPathsConfFileElement(includepath);
  else
    addParallelConfFileElement(inNode, inCluster, member_nodes, name,
                               computerName, port, minCpuNbr, maxCpuNbr, userName,
                               password, remoteDrive, remoteDirectory,
                               programPath, programConfig, matlabOctavePath, singleCompThread,
                               numberOfThreadsPerJob, operatingSystem);

  configFile.close();
}

void
ConfigFile::addHooksConfFileElement(string global_init_file)
{
  if (global_init_file.empty())
    {
      cerr << "ERROR: The global initialization file must be passed to the GlobalInitFile option." << endl;
      exit(EXIT_FAILURE);
    }
  else
    hooks.emplace_back(move(global_init_file));
}

void
ConfigFile::addPathsConfFileElement(vector<string> includepath)
{
  if (includepath.empty())
    {
      cerr << "ERROR: The path to be included must be passed to the Include option." << endl;
      exit(EXIT_FAILURE);
    }
  else
    paths.emplace_back(move(includepath));
}

void
ConfigFile::addParallelConfFileElement(bool inNode, bool inCluster, const member_nodes_t &member_nodes, const string &name,
                                       const string &computerName, const string &port, int minCpuNbr, int maxCpuNbr,
                                       const string &userName, const string &password, const string &remoteDrive,
                                       const string &remoteDirectory, const string &programPath, const string &programConfig,
                                       const string &matlabOctavePath, bool singleCompThread, int numberOfThreadsPerJob,
                                       const string &operatingSystem)
{
  //! ADD NODE
  if (inNode)
    if (!member_nodes.empty())
      {
        cerr << "Invalid option passed to [node]." << endl;
        exit(EXIT_FAILURE);
      }
    else
      if (name.empty() || follower_nodes.contains(name))
        {
          cerr << "ERROR: Every node must be assigned a unique name." << endl;
          exit(EXIT_FAILURE);
        }
      else
        follower_nodes.emplace(name, FollowerNode{computerName, port, minCpuNbr, maxCpuNbr, userName,
                                                  password, remoteDrive, remoteDirectory, programPath, programConfig,
                                                  matlabOctavePath, singleCompThread, numberOfThreadsPerJob,
                                                  operatingSystem});
  //! ADD CLUSTER
  else if (inCluster)
    if (minCpuNbr > 0 || maxCpuNbr > 0 || !userName.empty()
        || !password.empty() || !remoteDrive.empty() || !remoteDirectory.empty()
        || !programPath.empty() || !programConfig.empty()
        || !matlabOctavePath.empty() || !operatingSystem.empty())
      {
        cerr << "Invalid option passed to [cluster]." << endl;
        exit(EXIT_FAILURE);
      }
    else
      if (name.empty() || clusters.contains(name))
        {
          cerr << "ERROR: The cluster must be assigned a unique name." << endl;
          exit(EXIT_FAILURE);
        }
      else
        {
          if (clusters.empty())
            firstClusterName = name;
          clusters.emplace(name, Cluster{member_nodes});
        }
}

void
ConfigFile::checkPass(WarningConsolidation &warnings) const
{
  for (bool global_init_file_declared{false};
       const auto &hook : hooks)
    for (const auto &mapit : hook.get_hooks())
      if (mapit.first == "global_init_file")
        if (exchange(global_init_file_declared, true))
          {
            cerr << "ERROR: Only one global initialization file may be provided." << endl;
            exit(EXIT_FAILURE);
          }

  if (!parallel && !parallel_test)
    return;

  //! Check Follower Nodes
  if (follower_nodes.empty())
    {
      cerr << "ERROR: At least one node must be defined in the config file." << endl;
      exit(EXIT_FAILURE);
    }

  for (const auto &follower_node : follower_nodes)
    {
#if !defined(_WIN32) && !defined(__CYGWIN32__)
      //For Linux/Mac, check that cpuNbr starts at 0
      if (follower_node.second.minCpuNbr != 0)
        warnings << "WARNING: On Unix-based operating systems, you cannot specify the CPU that is "
                 << "used in parallel processing. This will be adjusted for you such that the "
                 << "same number of CPUs are used." << endl;
#endif
      if (!follower_node.second.port.empty())
        try
          {
            stoi(follower_node.second.port);
          }
        catch (const invalid_argument &)
          {
            cerr << "ERROR (node " << follower_node.first << "): the port must be an integer." << endl;
            exit(EXIT_FAILURE);
          }
      if (follower_node.second.computerName == "localhost") // We are working locally
        {
          if (!follower_node.second.remoteDrive.empty())
            {
              cerr << "ERROR (node " << follower_node.first << "): the RemoteDrive option may not be passed for a local node." << endl;
              exit(EXIT_FAILURE);
            }
          if (!follower_node.second.remoteDirectory.empty())
            {
              cerr << "ERROR (node " << follower_node.first << "): the RemoteDirectory option may not be passed for a local node." << endl;
              exit(EXIT_FAILURE);
            }
        }
      else
        {
          if (follower_node.second.userName.empty())
            {
              cerr << "ERROR (node " << follower_node.first << "): the UserName option must be passed for every remote node." << endl;
              exit(EXIT_FAILURE);
            }
          if (follower_node.second.operatingSystem == "windows")
            {
              if (follower_node.second.password.empty())
                {
                  cerr << "ERROR (node " << follower_node.first << "): the Password option must be passed under Windows for every remote node." << endl;
                  exit(EXIT_FAILURE);
                }
              if (follower_node.second.remoteDrive.empty())
                {
                  cerr << "ERROR (node " << follower_node.first << "): the RemoteDrive option must be passed under Windows for every remote node." << endl;
                  exit(EXIT_FAILURE);
                }
            }
#if defined(_WIN32) || defined(__CYGWIN32__)
          if (follower_node.second.operatingSystem.empty())
            {
              if (follower_node.second.password.empty())
                {
                  cerr << "ERROR (node " << follower_node.first << "): the Password option must be passed under Windows for every remote node." << endl;
                  exit(EXIT_FAILURE);
                }
              if (follower_node.second.remoteDrive.empty())
                {
                  cerr << "ERROR (node " << follower_node.first << "): the RemoteDrive option must be passed under Windows for every remote node." << endl;
                  exit(EXIT_FAILURE);
                }
            }
#endif
          if (follower_node.second.remoteDirectory.empty())
            {
              cerr << "ERROR (node " << follower_node.first << "): the RemoteDirectory must be specified for every remote node." << endl;
              exit(EXIT_FAILURE);
            }
        }
    }

  //! Check Clusters
  if (clusters.empty())
    {
      cerr << "ERROR: At least one cluster must be defined in the config file." << endl;
      exit(EXIT_FAILURE);
    }

  if (!cluster_name.empty() && !clusters.contains(cluster_name))
    {
      cerr << "ERROR: Cluster Name " << cluster_name << " was not found in the config file." << endl;
      exit(EXIT_FAILURE);
    }

  for (const auto &cluster : clusters)
    for (const auto &itmn : cluster.second.member_nodes)
      if (!follower_nodes.contains(itmn.first))
        {
          cerr << "Error: node " << itmn.first << " specified in cluster " << cluster.first << " was not found" << endl;
          exit(EXIT_FAILURE);
        }
}

void
ConfigFile::transformPass()
{
  if (!parallel && !parallel_test)
    return;

#if !defined(_WIN32) && !defined(__CYGWIN32__)
  //For Linux/Mac, check that cpuNbr starts at 0
  for (auto &it : follower_nodes)
    if (it.second.minCpuNbr != 0)
      {
        it.second.maxCpuNbr = it.second.maxCpuNbr - it.second.minCpuNbr;
        it.second.minCpuNbr = 0;
      }
#endif

  auto cluster_it = cluster_name.empty() ? clusters.find(firstClusterName) : clusters.find(cluster_name);

  double weight_denominator{0.0};
  for (const auto &it : cluster_it->second.member_nodes)
    weight_denominator += it.second;

  for (auto &member_node : cluster_it->second.member_nodes)
    member_node.second /= weight_denominator;
}

vector<filesystem::path>
ConfigFile::getIncludePaths() const
{
  vector<filesystem::path> include_paths;
  for (auto path : paths)
    for (const auto &mapit : path.get_paths())
      for (const auto &vecit : mapit.second)
        include_paths.emplace_back(vecit);
  return include_paths;
}

void
ConfigFile::writeHooks(ostream &output) const
{
  for (auto hook : hooks)
    for (const auto &mapit : hook.get_hooks())
      output << "options_." << mapit.first << " = '" << mapit.second << "';" << endl;
}

void
ConfigFile::writeCluster(ostream &output) const
{
  if (!parallel && !parallel_test)
    return;

  auto cluster_it = cluster_name.empty() ? clusters.find(firstClusterName) : clusters.find(cluster_name);

  for (int i{1};
       const auto &follower_node : follower_nodes)
    {
      bool follower_node_in_member_nodes = false;
      for (const auto &itmn : cluster_it->second.member_nodes)
        if (follower_node.first == itmn.first)
          follower_node_in_member_nodes = true;

      if (!follower_node_in_member_nodes)
        continue;

      output << "options_.parallel";
      if (i > 1)
        output << "(" << i << ")";
      i++;
      output << " = struct('Local', ";
      if (follower_node.second.computerName == "localhost")
        output << "1, ";
      else
        output << "0, ";

      output << "'ComputerName', '" << follower_node.second.computerName << "', "
             << "'Port', '" << follower_node.second.port << "', "
             << "'CPUnbr', [" << follower_node.second.minCpuNbr << ":" << follower_node.second.maxCpuNbr << "], "
             << "'UserName', '" << follower_node.second.userName << "', "
             << "'Password', '" << follower_node.second.password << "', "
             << "'RemoteDrive', '" << follower_node.second.remoteDrive << "', "
             << "'RemoteDirectory', '" << follower_node.second.remoteDirectory << "', "
        // The following should be switched back to “ProgramPath” once we move to Dragonfly
             << "'DynarePath', '" << follower_node.second.programPath << "', "
             << "'ProgramConfig', '" << follower_node.second.programConfig << "', "
             << "'MatlabOctavePath', '" << follower_node.second.matlabOctavePath << "', "
             << "'OperatingSystem', '" << follower_node.second.operatingSystem << "', "
             << "'NodeWeight', '" << cluster_it->second.member_nodes.at(follower_node.first) << "', "
             << "'NumberOfThreadsPerJob', " << follower_node.second.numberOfThreadsPerJob << ", ";

      if (follower_node.second.singleCompThread)
        output << "'SingleCompThread', 'true');" << endl;
      else
        output << "'SingleCompThread', 'false');" << endl;
    }

  // Default values for the following two are both in DynareMain.cc and matlab/default_option_values.m
  if (parallel_follower_open_mode)
    output << "options_.parallel_info.leaveSlaveOpen = 1;" << endl;
  if (!parallel_use_psexec)
    output << "options_.parallel_use_psexec = false;" << endl;

  output << "options_.parallel_info.console_mode= isoctave;" << endl;

  output << "InitializeComputationalEnvironment();" << endl;
  if (parallel_test)
    output << "ErrorCode = AnalyseComputationalEnvironment(options_.parallel, options_.parallel_info);" << endl
           << "disp(['AnalyseComputationalEnvironment returned with Error Code: ' num2str(ErrorCode)]);" << endl
           << "diary off;" << endl
           << "return;" << endl;
}

void
ConfigFile::writeEndParallel(ostream &output) const
{
  if ((!parallel && !parallel_test) || !parallel_follower_open_mode)
    return;

  output << "if options_.parallel_info.leaveSlaveOpen == 1" << endl
         << "     closeSlave(options_.parallel,options_.parallel_info.RemoteTmpFolder);" << endl
         << "end" << endl;
}
