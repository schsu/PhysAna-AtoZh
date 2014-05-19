#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <cstring>

#include <unistd.h>
#include <sys/wait.h>
#include <errno.h>
std::map<int, std::string> runs;

void run_pythia();
void run_pythia(std::string mass);

char* strdup(const char* s) {
  char *p = (char*)malloc(std::strlen(s) + 1);
  std::strcpy(p, s);
  return p;
}

void replace_template(std::string& line, const std::string& from, const std::string& to) {
  size_t am = line.find(from);
  if (am != std::string::npos) {
    line.replace(am, from.length(), to);
  }
}

void expand_range(int from, int to, const std::string& input, const std::string& output) {
  std::ifstream inp(input.c_str());
  std::ofstream outp(output.c_str(), std::ios::out);
	
  std::cout << input.c_str() << " -> " << output.c_str() << std::endl;

  char buffer[32];
  sprintf(buffer, "%d", from);
  std::string froms(buffer);
  sprintf(buffer, "%d", to);
  std::string tos(buffer);

  while (!inp.eof()) {
    std::string line;
    std::getline(inp, line);
    replace_template(line, "$ihtmin$", froms);
    replace_template(line, "$ihtmax$", tos);
    outp << line << std::endl;
  }
}

void expand_mass(const std::string& mass, const std::string& input, const std::string& output) {
  std::ifstream inp(input.c_str());
  std::ofstream outp(output.c_str(), std::ios::out);

  while (!inp.eof()) {
    std::string line;
    std::getline(inp, line);
    replace_template(line, "$amass$", mass);
    outp << line << std::endl;
  }
}

void execute_cmd_no_wait(const std::string& cmd, std::vector<std::string>& params) {
  pid_t procId = fork();
  if (procId == 0) {
    char** args = new char*[params.size() + 2];
    args[0] = strdup(cmd.c_str());
    args[params.size()+1] = NULL;
    for (size_t i = 0; i < params.size(); ++i) {
      args[i+1] = strdup(params[i].c_str());
    }

    execvp(cmd.c_str(), args);
  }
}

void execute_cmd_wait(const std::string& cmd, std::vector<std::string>& params) {
  pid_t procId = fork();
  // procId == 0 -> we forked, this is the new process
  if (procId == 0) {
    char** args = new char*[params.size() + 2];
    args[0] = strdup(cmd.c_str());
    args[params.size()+1] = NULL;
    for (size_t i = 0; i < params.size(); ++i) {
      args[i+1] = strdup(params[i].c_str());
    }

    execvp(cmd.c_str(), args);
    std::cout << "failed. " << errno << std::endl;
  } else if (procId > 0) {
    int status = 0;
    do {
      waitpid(procId, &status, 0);
    } while (!WIFEXITED(status));
  } else {
    std::cerr << "Failed to create process." << std::endl;
  }
}

std::string get_mdg_root() {
  static std::string mdg;
  if (!mdg.empty()) {
    return mdg;
  }

  char* mdg_ch = std::getenv("MDGROOT");
  if (mdg_ch == NULL) {
    std::cerr << "Need to set MDGROOT." << std::endl;
    return "";
  }

  mdg = mdg_ch;
  return mdg;
}

void run_delphes(int higgsWeight, const std::string& procDirectory) {
  char buffer[32];
  std::vector<std::string> params;

  std::string delphesPath = get_mdg_root() + "Delphes/DelphesSTDHEP";
  std::string cardPath = get_mdg_root() + "proc-2HDM4TC-gg-h3-zh/Cards/delphes-card-dual.dat";
  sprintf(buffer, "%d", higgsWeight);
  std::string outputPath = std::string("/phys/groups/tev/scratch3/users/milenl/A-Zh/outputs/a-zh-") + buffer + "GeV.root";
  std::string inputPath = get_mdg_root() + procDirectory + "/Events/" + runs[higgsWeight] + "/tag_1_pythia_events.hep";

  std::remove(outputPath.c_str());

  params.push_back(cardPath);
  params.push_back(outputPath);
  params.push_back(inputPath);

  std::cout << delphesPath.c_str() << std::endl;
  execute_cmd_wait(delphesPath, params);
}

void initialize_heavy_h() {
  runs[250] = "run_01";
  runs[300] = "run_02";
  runs[350] = "run_03";
  runs[400] = "run_04";
  runs[500] = "run_05";
  runs[600] = "run_06";
  runs[700] = "run_07";
  runs[800] = "run_08";
  runs[900] = "run_09";
  runs[1000] = "run_10";
}

void run_pythia() {
  run_pythia("2.5");
  run_pythia("3.0");
  run_pythia("3.5");
  run_pythia("4.0");
  run_pythia("5.0");
  run_pythia("6.0");
  run_pythia("7.0");
  run_pythia("8.0");
  run_pythia("9.0");
  run_pythia("10.0");
  run_pythia("11.0");
  run_pythia("12.0");
  run_pythia("13.0");
  run_pythia("14.0");
  run_pythia("15.0");
}

void initialize_background() {
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 8; ++j) {
      char buff[2048];
      sprintf(buff, "run_%02d", i*9+j+1);
      runs[101+i*9+j] = strdup(buff);
    }
  }
}

void run_pythia(std::string mass) {
  std::string paramTemplate = get_mdg_root() + "proc-2HDM4TC-gg-h3-zh/Cards/param_card_template.dat";
  std::string paramCard = get_mdg_root() + "proc-2HDM4TC-gg-h3-zh/Cards/param_card.dat";
  expand_mass(mass, paramTemplate, paramCard);

  std::string pythiaPath = get_mdg_root() + "proc-2HDM4TC-gg-h3-zh/bin/generate_events";
  std::vector<std::string> params;
  params.push_back("--laststep=pythia");
  execute_cmd_wait(pythiaPath, params);
}

void run_background_me(int from, int to) {
  std::string paramTemplate = get_mdg_root() + "proc-std-pp-zbb/Cards/run_card_template.dat";
  std::string paramRunCard = get_mdg_root() + "proc-std-pp-zbb/Cards/run_card.dat";
  expand_range(from, to, paramTemplate, paramRunCard);
  std::string madEventPath = get_mdg_root() + "proc-std-pp-zbb/bin/generate_events";
  std::vector<std::string> params;
  params.push_back("--laststep=pythia");
  execute_cmd_wait(madEventPath, params);
}

void run_background_me() {
  run_background_me(0, 100);
  run_background_me(100, 200);
  run_background_me(200, 300);
  run_background_me(300, 400);
  run_background_me(400, 500);
  run_background_me(500, 600);
  run_background_me(600, 800);
  run_background_me(800, 1000);
  run_background_me(300,800);
}

void
generate_events_signal() {
  run_pythia();
}

void remote_execute_cmd(const std::string& machine, const std::string& cmd) {
  std::vector<std::string> params;
  params.push_back(machine);
  params.push_back(cmd);
  execute_cmd_wait("ssh", params);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: batchrun -generate-signal|-generate-background|-delphes-signal|-delphes-background" << std::endl;
    return 2;
  }
  
  char* mdg_ch = std::getenv("MDGROOT");
  if (mdg_ch == NULL) {
    std::cerr << "Need to set MDGROOT. (like export MDGROOT=/path/to/MadGraph5/)" << std::endl;
    return 1;
  }

  if (strcmp(argv[1], "-generate-signal") == 0) {
    run_pythia();
  } else if (strcmp(argv[1], "-generate-background") == 0) {
    run_background_me();
  } else if (strcmp(argv[1], "-delphes-signal") == 0) {
    initialize_heavy_h();
    for (std::map<int, std::string>::iterator i = runs.begin(); i != runs.end(); ++i) {
      run_delphes(i->first, "proc-2HDM4TC-gg-h3-zh");
    }
  } else if (strcmp(argv[1], "-delphes-background") == 0) {
    initialize_background();
    for (std::map<int, std::string>::iterator i = runs.begin(); i != runs.end(); ++i) {
      run_delphes(i->first, "proc-std-pp-zbb");
    }
  } else {
    std::cout << "Usage: batchrun -generate-signal|-generate-background|-delphes-signal|-delphes-background" << std::endl;
  }
 
  return 0;
}
