/*
 * Copyright © 2010-2022 Dynare Team
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

#ifndef _CONFIG_FILE_HH
#define _CONFIG_FILE_HH

#include <map>
#include <vector>

#include "WarningConsolidation.hh"

using namespace std;

using member_nodes_t = map<string, double>;

class Hook
{
public:
  explicit Hook(string global_init_file_arg);
private:
  map<string, string> hooks;
public:
  map<string, string>
  get_hooks() const
  {
    return hooks;
  };
};

class Path
{
public:
  explicit Path(vector<string> includepath_arg);
private:
  map<string, vector<string>> paths;
public:
  map<string, vector<string>>
  get_paths() const
  {
    return paths;
  };
};

class FollowerNode
{
  friend class ConfigFile;
public:
  FollowerNode(string computerName_arg, string port_arg, int minCpuNbr_arg, int maxCpuNbr_arg, string userName_arg,
               string password_arg, string remoteDrive_arg, string remoteDirectory_arg,
               string programPath_arg, string programConfig_arg, string matlabOctavePath_arg, bool singleCompThread_arg,
               int numberOfThreadsPerJob_arg, string operatingSystem_arg);

protected:
  const string computerName, port;
  int minCpuNbr, maxCpuNbr;
  const string userName, password;
  const string remoteDrive, remoteDirectory;
  const string programPath, programConfig, matlabOctavePath;
  const bool singleCompThread;
  const int numberOfThreadsPerJob;
  const string operatingSystem;
};

class Cluster
{
  friend class ConfigFile;
public:
  explicit Cluster(member_nodes_t member_nodes_arg);

protected:
  member_nodes_t member_nodes;
};

//! The abstract representation of a "config" file
class ConfigFile
{
public:
  ConfigFile(bool parallel_arg, bool parallel_test_arg, bool parallel_follower_open_mode_arg,
             bool parallel_use_psexec_arg, string cluster_name);

private:
  const bool parallel, parallel_test, parallel_follower_open_mode, parallel_use_psexec;
  const string cluster_name;
  string firstClusterName;
  //! Hooks
  vector<Hook> hooks;
  //! Paths
  vector<Path> paths;
  //! Cluster Table
  map<string, Cluster> clusters;
  //! Node Map
  map<string, FollowerNode> follower_nodes;
  //! Add Hooks
  void addHooksConfFileElement(string global_init_file);
  //! Add Paths
  void addPathsConfFileElement(vector<string> includepath);
  //! Add a FollowerNode or a Cluster object
  void addParallelConfFileElement(bool inNode, bool inCluster, const member_nodes_t &member_nodes, const string &name,
                                  const string &computerName, const string &port, int minCpuNbr, int maxCpuNbr,
                                  const string &userName, const string &password, const string &remoteDrive,
                                  const string &remoteDirectory, const string &programPath, const string &programConfig,
                                  const string &matlabOctavePath, bool singleCompThread, int numberOfThreadsPerJob,
                                  const string &operatingSystem);
public:
  //! Parse config file
  void getConfigFileInfo(const string &parallel_config_file);
  //! Check Pass
  void checkPass(WarningConsolidation &warnings) const;
  //! Check Pass
  void transformPass();
  //! Get Path Info
  vector<filesystem::path> getIncludePaths() const;
  //! Write any hooks
  void writeHooks(ostream &output) const;
  //! Create options_.parallel structure, write options
  void writeCluster(ostream &output) const;
  //! Close follower nodes if needed
  void writeEndParallel(ostream &output) const;
};

#endif // ! CONFIG_FILE_HH
