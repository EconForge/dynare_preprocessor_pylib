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

#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <filesystem>
#include <map>
#include <vector>

#include "WarningConsolidation.hh"

using namespace std;

/* The abstract representation of the configuration.
   Merges information from the command-line and from the configuration file. */
class Configuration
{
public:
  Configuration(bool parallel_arg, bool parallel_test_arg, bool parallel_follower_open_mode_arg,
                bool parallel_use_psexec_arg, string cluster_name);

private:
  using member_nodes_t = map<string, double>;

  class Path
  {
  public:
    explicit Path(vector<string> includepath_arg);
    [[nodiscard]] map<string, vector<string>>
    get_paths() const
    {
      return paths;
    };

  private:
    map<string, vector<string>> paths;
  };

  struct FollowerNode
  {
    FollowerNode(string computerName_arg, string port_arg, int minCpuNbr_arg, int maxCpuNbr_arg,
                 string userName_arg, string password_arg, string remoteDrive_arg,
                 string remoteDirectory_arg, string programPath_arg, string programConfig_arg,
                 string matlabOctavePath_arg, bool singleCompThread_arg,
                 int numberOfThreadsPerJob_arg, string operatingSystem_arg);
    const string computerName, port;
    int minCpuNbr, maxCpuNbr;
    const string userName, password;
    const string remoteDrive, remoteDirectory;
    const string programPath, programConfig, matlabOctavePath;
    const bool singleCompThread;
    const int numberOfThreadsPerJob;
    const string operatingSystem;
  };

  struct Cluster
  {
    explicit Cluster(member_nodes_t member_nodes_arg);
    member_nodes_t member_nodes;
  };

  const bool parallel, parallel_test, parallel_follower_open_mode, parallel_use_psexec;
  const string cluster_name;
  string firstClusterName;
  //! Hooks
  string global_init_file;
  //! Paths
  vector<Path> paths;
  //! Cluster Table
  map<string, Cluster> clusters;
  //! Node Map
  map<string, FollowerNode> follower_nodes;
  //! Add Paths
  void addPathsConfFileElement(vector<string> includepath);
  //! Add a FollowerNode or a Cluster object
  void addParallelConfFileElement(bool inNode, bool inCluster, const member_nodes_t& member_nodes,
                                  const string& name, const string& computerName,
                                  const string& port, int minCpuNbr, int maxCpuNbr,
                                  const string& userName, const string& password,
                                  const string& remoteDrive, const string& remoteDirectory,
                                  const string& programPath, const string& programConfig,
                                  const string& matlabOctavePath, bool singleCompThread,
                                  int numberOfThreadsPerJob, const string& operatingSystem);
  /* Given a filename (e.g. dynare.ini), looks for it in the configuration directory:
     – if under Linux or macOS, look into the “dynare” subdirectory of the XDG
       configuration directories (following the default values and the precedence order specified in
       the XDG specification)
     – if under Windows, look into %APPDATA%\dynare\
     The returned path will be empty if the file is not found. */
  [[nodiscard]] static filesystem::path findConfigFile(const string& filename);

public:
  //! Parse config file
  void getConfigFileInfo(const filesystem::path& conffile_option, WarningConsolidation& warnings);
  //! Check Pass
  void checkPass(WarningConsolidation& warnings) const;
  //! Check Pass
  void transformPass();
  //! Get Path Info
  [[nodiscard]] vector<filesystem::path> getIncludePaths() const;
  //! Write any hooks
  void writeHooks(ostream& output) const;
  //! Create options_.parallel structure, write options
  void writeCluster(ostream& output) const;
  //! Close follower nodes if needed
  void writeEndParallel(ostream& output) const;
};

#endif
