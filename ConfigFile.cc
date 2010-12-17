/*
 * Copyright (C) 2010 Dynare Team
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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "ConfigFile.hh"
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;

SlaveNode::SlaveNode(string &computerName_arg, int minCpuNbr_arg, int maxCpuNbr_arg, string &userName_arg,
                     string &password_arg, string &remoteDrive_arg, string &remoteDirectory_arg,
                     string &dynarePath_arg,  string &matlabOctavePath_arg, bool singleCompThread_arg) :
  computerName(computerName_arg), minCpuNbr(minCpuNbr_arg), maxCpuNbr(maxCpuNbr_arg), userName(userName_arg),
  password(password_arg), remoteDrive(remoteDrive_arg), remoteDirectory(remoteDirectory_arg), dynarePath(dynarePath_arg),
  matlabOctavePath(matlabOctavePath_arg), singleCompThread(singleCompThread_arg)
{
  if (computerName.empty())
    {
      cerr << "ERROR: The node must have a ComputerName." << endl;
      exit(EXIT_FAILURE);
    }
}

Cluster::Cluster(vector<string> member_nodes_arg) : member_nodes(member_nodes_arg)
{
  if (member_nodes.empty())
    {
      cerr << "ERROR: The cluster must have at least one member node." << endl;
      exit(EXIT_FAILURE);
    }
}

ConfigFile::ConfigFile(bool parallel_arg, bool parallel_test_arg,
                       bool parallel_slave_open_mode_arg, const string &cluster_name_arg) :
  parallel(parallel_arg), parallel_test(parallel_test_arg),
  parallel_slave_open_mode(parallel_slave_open_mode_arg), cluster_name(cluster_name_arg)
{
}

ConfigFile::~ConfigFile()
{
}

void
ConfigFile::getConfigFileInfo(const string &parallel_config_file)
{
  if (!parallel && !parallel_test)
    return;

  ifstream *configFile;
  if (parallel_config_file.empty())
    {
      // Test OS and try to open default file
#if defined(_WIN32) || defined(__CYGWIN32__)
      if (getenv("APPDATA")==NULL)
        {
          cerr << "ERROR: APPDATA environment variable not found." << endl;
          exit(EXIT_FAILURE);
        }
      string defaultConfigFile (getenv("APPDATA"));
      defaultConfigFile += "\\dynare.ini";
#else
      if (getenv("HOME")==NULL)
        {
          cerr << "ERROR: HOME environment variable not found." << endl;
          exit(EXIT_FAILURE);
        }
      string defaultConfigFile (getenv("HOME"));
      defaultConfigFile += "/.dynare";
#endif
      configFile = new ifstream(defaultConfigFile.c_str(), fstream::in);
      if (!configFile->is_open())
        {
          cerr << "ERROR: Could not open the default config file (" << defaultConfigFile << ")" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      configFile = new ifstream(parallel_config_file.c_str(), fstream::in);
      if (!configFile->is_open())
        {
          cerr << "ERROR: Couldn't open file " << parallel_config_file << endl;;
          exit(EXIT_FAILURE);
        }
    }

  string name, computerName, userName, password, remoteDrive,
    remoteDirectory, dynarePath, matlabOctavePath;
  int minCpuNbr = 0, maxCpuNbr = 0;
  bool singleCompThread = true;
  vector<string> member_nodes;

  bool inNode = false;
  bool inCluster = false;
  while (configFile->good())
    {
      string line;
      getline(*configFile, line);
      boost::trim(line);
      if (line.empty())
        continue;

      if (!line.compare("[node]") || !line.compare("[cluster]"))
        {
          addConfFileElement(inNode, inCluster, member_nodes, name,
                 computerName, minCpuNbr, maxCpuNbr, userName,
                 password, remoteDrive, remoteDirectory,
                 dynarePath,  matlabOctavePath, singleCompThread);

          //! Reset communication vars / option defaults
          if (!line.compare("[node]"))
            {
              inNode = true;
              inCluster = false;
            }
          else
            {
              inNode = false;
              inCluster = true;
            }

          name = userName = computerName = password = remoteDrive =
            remoteDirectory = dynarePath = matlabOctavePath = "";
          minCpuNbr = maxCpuNbr = 0;
          singleCompThread = true;
          member_nodes.clear();
        }
      else
        {
          vector<string> tokenizedLine;
          boost::split(tokenizedLine, line, boost::is_any_of("="));
          if (tokenizedLine.size() != 2)
            {
              cerr << "ERROR (in config file): Options should be formatted as 'option = value'." << endl;
              exit(EXIT_FAILURE);
            }
          boost::trim(tokenizedLine.front());
          boost::trim(tokenizedLine.back());

          if (!tokenizedLine.front().compare("Name"))
            name = tokenizedLine.back();
          else if (!tokenizedLine.front().compare("CPUnbr"))
            {
              vector<string> tokenizedCpuNbr;
              boost::split(tokenizedCpuNbr, tokenizedLine.back(), boost::is_any_of(":"));
              try
                {
                  if (tokenizedCpuNbr.size() == 1)
                    {
                      minCpuNbr = 1;
                      maxCpuNbr = boost::lexical_cast< int >(tokenizedCpuNbr.front());
                    }
                  else if (tokenizedCpuNbr.size() == 2 &&
                           tokenizedCpuNbr[0].at(0) == '[' &&
                           tokenizedCpuNbr[1].at(tokenizedCpuNbr[1].size()-1) == ']')
                    {
                      tokenizedCpuNbr[0].erase(0,1);
                      tokenizedCpuNbr[1].erase(tokenizedCpuNbr[1].size()-1,1);
                      minCpuNbr = boost::lexical_cast< int >(tokenizedCpuNbr[0]);
                      maxCpuNbr = boost::lexical_cast< int >(tokenizedCpuNbr[1]);
                    }
                }
              catch( const boost::bad_lexical_cast & )
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
          else if (!tokenizedLine.front().compare("ComputerName"))
            computerName = tokenizedLine.back();
          else if (!tokenizedLine.front().compare("UserName"))
            userName = tokenizedLine.back();
          else if (!tokenizedLine.front().compare("Password"))
            password = tokenizedLine.back();
          else if (!tokenizedLine.front().compare("RemoteDrive"))
            remoteDrive = tokenizedLine.back();
          else if (!tokenizedLine.front().compare("RemoteDirectory"))
            remoteDirectory = tokenizedLine.back();
          else if (!tokenizedLine.front().compare("DynarePath"))
            dynarePath = tokenizedLine.back();
          else if (!tokenizedLine.front().compare("MatlabOctavePath"))
            matlabOctavePath = tokenizedLine.back();
          else if (!tokenizedLine.front().compare("SingleCompThread"))
            if (tokenizedLine.back().compare("true") == 0)
              singleCompThread = true;
            else if (tokenizedLine.back().compare("false") == 0)
              singleCompThread = false;
            else
              {
                cerr << "ERROR (in config file): The value passed to SingleCompThread may only be 'true' or 'false'." << endl;
                exit(EXIT_FAILURE);
              }
          else if (!tokenizedLine.front().compare("Members"))
            {
              vector<string> tmp_member_nodes;
              boost::split(tmp_member_nodes, tokenizedLine.back(), boost::is_any_of(";, "));
              for ( vector<string>::iterator it = tmp_member_nodes.begin();
                    it < tmp_member_nodes.end(); it++ )
                {
                  boost::trim(*it);
                  if (!it->empty())
                    member_nodes.push_back(*it);
                }
            }
          else
            {
              cerr << "ERROR (in config file): Option " << tokenizedLine.front() << " is invalid." << endl;
              exit(EXIT_FAILURE);
            }
        }
    }

  addConfFileElement(inNode, inCluster, member_nodes, name,
                     computerName, minCpuNbr, maxCpuNbr, userName,
                     password, remoteDrive, remoteDirectory,
                     dynarePath,  matlabOctavePath, singleCompThread);
  configFile->close();
  delete configFile;
}

void
ConfigFile::addConfFileElement(bool inNode, bool inCluster, vector<string> member_nodes, string &name,
                               string &computerName, int minCpuNbr, int maxCpuNbr, string &userName,
                               string &password, string &remoteDrive, string &remoteDirectory,
                               string &dynarePath,  string &matlabOctavePath, bool singleCompThread)
{
  //! ADD NODE
  if (inNode)
    if (!member_nodes.empty())
      {
        cerr << "Invalid option passed to [node]." << endl;
        exit(EXIT_FAILURE);
      }
    else
      if (name.empty() || slave_nodes.find(name) != slave_nodes.end())
        {
          cerr << "ERROR: Every node must be assigned a unique name." << endl;
          exit(EXIT_FAILURE);
        }
      else
        slave_nodes[name] = new SlaveNode(computerName, minCpuNbr, maxCpuNbr, userName,
                                          password, remoteDrive, remoteDirectory, dynarePath,
                                          matlabOctavePath, singleCompThread);
  //! ADD CLUSTER
  else if (inCluster)
    if ( minCpuNbr > 0 || maxCpuNbr > 0 || !userName.empty() ||
         !password.empty() || !remoteDrive.empty() || !remoteDirectory.empty() ||
         !dynarePath.empty() || !matlabOctavePath.empty())
      {
        cerr << "Invalid option passed to [cluster]." << endl;
        exit(EXIT_FAILURE);
      }
    else
      if (name.empty() || clusters.find(name) != clusters.end())
        {
          cerr << "ERROR: The cluster must be assigned a unique name." << endl;
          exit(EXIT_FAILURE);
        }
      else
        {
          if (clusters.empty())
            firstClusterName = name;
          clusters[name] = new Cluster(member_nodes);
        }
}

void
ConfigFile::checkPass() const
{
  if (!parallel && !parallel_test)
    return;

  //! Check Slave Nodes
  if (slave_nodes.empty())
    {
      cerr << "ERROR: At least one node must be defined in the config file." << endl;
      exit(EXIT_FAILURE);
    }

  for (map<string, SlaveNode *>::const_iterator it = slave_nodes.begin();
       it != slave_nodes.end(); it++)
    {
#if !defined(_WIN32) && !defined(__CYGWIN32__)
      //For Linux/Mac, check that cpuNbr starts at 0
      if (it->second->minCpuNbr != 0)
        cout << "WARNING: On Unix-based operating systems, you cannot specify the CPU that is used "
             << "in parallel processing. This will be adjusted for you such that the same number of CPUs "
             << "are used." << endl;
#endif
      if (!it->second->computerName.compare("localhost")) // We are working locally
        {
          if (!it->second->remoteDrive.empty())
            {
              cerr << "ERROR (node " << it->first << "): the RemoteDrive option may not be passed for a local node." << endl;
              exit(EXIT_FAILURE);
            }
          if (!it->second->remoteDirectory.empty())
            {
              cerr << "ERROR (node " << it->first << "): the RemoteDirectory option may not be passed for a local node." << endl;
              exit(EXIT_FAILURE);
            }
        }
      else
        {
          if (it->second->userName.empty())
            {
              cerr << "ERROR (node " << it->first << "): the UserName option must be passed for every remote node." << endl;
              exit(EXIT_FAILURE);
            }
#if defined(_WIN32) || defined(__CYGWIN32__)
          if (it->second->userName.empty() || it->second->password.empty())
            {
              cerr << "ERROR (node " << it->first << "): the Password option must be passed under Windows for every remote node." << endl;
              exit(EXIT_FAILURE);
            }
          if (it->second->remoteDrive.empty())
            {
              cerr << "ERROR (node " << it->first << "): the RemoteDrive option must be passed under Windows for every remote node." << endl;
              exit(EXIT_FAILURE);
            }
#endif
          if (it->second->remoteDirectory.empty())
            {
              cerr << "ERROR (node " << it->first << "): the RemoteDirectory must be specified for every remote node." << endl;
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

  if (!cluster_name.empty() && clusters.find(cluster_name) == clusters.end())
    {
      cerr << "ERROR: Cluster Name " << cluster_name << " was not found in the config file." << endl;
      exit(EXIT_FAILURE);
    }

  for (map<string, Cluster *>::const_iterator it = clusters.begin();
       it != clusters.end(); it++)
    for (vector<string>::const_iterator itmn = it->second->member_nodes.begin();
         itmn < it->second->member_nodes.end(); itmn++)
      if (slave_nodes.find(*itmn) == slave_nodes.end())
        {
          cerr << "Error: node " << *itmn << " specified in cluster " << it->first << " was not found" << endl;
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
  for (map<string, SlaveNode *>::const_iterator it = slave_nodes.begin();
       it != slave_nodes.end(); it++)
    if (it->second->minCpuNbr != 0)
      {
        it->second->maxCpuNbr = it->second->maxCpuNbr - it->second->minCpuNbr;
        it->second->minCpuNbr = 0;
      }
#endif
}

void
ConfigFile::writeCluster(ostream &output) const
{
  if (!parallel && !parallel_test)
    return;

  map<string, Cluster *>::const_iterator cluster_it ;
  if (cluster_name.empty())
    cluster_it = clusters.find(firstClusterName);
  else
    cluster_it = clusters.find(cluster_name);

  int i = 1;
  for (map<string, SlaveNode *>::const_iterator it = slave_nodes.begin();
       it != slave_nodes.end(); it++)
    {
      bool slave_node_in_member_nodes = false;
      for (vector<string>::const_iterator itmn = cluster_it->second->member_nodes.begin();
           itmn < cluster_it->second->member_nodes.end(); itmn++)
        if (!it->first.compare(*itmn))
          slave_node_in_member_nodes = true;

      if (!slave_node_in_member_nodes)
        continue;

      output << "options_.parallel";
      if (i > 1)
        output << "(" << i << ")";
      i++;
      output << " = struct('Local', ";
      if (it->second->computerName.compare("localhost"))
        output << "0, ";
      else
        output << "1, ";

      output << "'ComputerName', '" << it->second->computerName << "', "
             << "'CPUnbr', [" << it->second->minCpuNbr << ":" << it->second->maxCpuNbr << "], "
             << "'UserName', '" << it->second->userName << "', "
             << "'Password', '" << it->second->password << "', "
             << "'RemoteDrive', '" << it->second->remoteDrive << "', "
             << "'RemoteDirectory', '" << it->second->remoteDirectory << "', "
             << "'DynarePath', '" << it->second->dynarePath << "', "
             << "'MatlabOctavePath', '" << it->second->matlabOctavePath << "', ";

      if (it->second->singleCompThread)
        output << "'SingleCompThread', 'true');" << endl;
      else
        output << "'SingleCompThread', 'false');" << endl;
    }

  if (parallel_slave_open_mode)
    output << "options_.parallel_info.leaveSlaveOpen = 1;" << endl;

  output << "InitializeComputationalEnvironment();" << endl;
  if (parallel_test)
    output  << "ErrorCode = AnalyseComputationalEnvironment(options_.parallel, options_.parallel_info);" << endl
            << "disp(['AnalyseComputationalEnvironment returned with Error Code: ' num2str(ErrorCode)]);" << endl
            << "diary off;" << endl
            << "return;" << endl;
}

void
ConfigFile::writeEndParallel(ostream &output) const
{
  if ((!parallel && !parallel_test) || !parallel_slave_open_mode)
    return;

  output << "if options_.parallel_info.leaveSlaveOpen == 1" << endl
         << "     closeSlave(options_.parallel,options_.parallel_info.RemoteTmpFolder);" << endl
         << "end" << endl;
}
