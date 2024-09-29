// Extract data from clinvar.vcf based on rsIDs found in AncestryDNA.txt
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <regex>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// Structure to hold ClinVar data
struct ClinVarData {
    std::string chrom;
    std::string pos;
    std::string id;
    std::string refAllele;
    std::string altAllele;
    std::string qual;
    std::string filter;
    std::string geneInfo;
    std::string info;
};

// Helper function to split a string by a delimiter
std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> result;
    std::istringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        result.push_back(token);
    }
    return result;
}
void trim_endlines(std::string &str) {
    // Keep removing characters while the last character is \r or \n
    while (!str.empty() && (str.back() == '\r' || str.back() == '\n')) {
        str.pop_back();  // Remove the last character
    }
}

// Function to process list of , separated items
json processKeyValue2( const std::string& value) {
    json jsonList;
    json jsonDict;
    // Find values or KEY:value pairs separated by ,
    std::vector<std::string> items = splitString(value, ',');
    for (const auto& item : items) {
        std::vector<std::string> kv = splitString(item, ':');
        if (kv.size() == 2) {
          jsonDict = json();
          jsonDict[kv[0]] = kv[1];
          jsonList += jsonDict;
        } else {
          jsonList += kv[0];
        }
    }
    return jsonList;
}

// Function to process list of | separated items
json processKeyValue1( const std::string& value) {
    json jsonData;
    // Find KEY=value pairs separated by ;
    std::vector<std::string> items = splitString(value, '|');
    for (const auto& item : items) {
      jsonData += processKeyValue2(item);
    }
    return jsonData;
}

// Recursive function to process key=value pairs and add them to a JSON object
void processKeyValue0(json& jsonData, const std::string& inputLine) {
    // Find KEY=value pairs separated by ;
    std::vector<std::string> keyValuePairs = splitString(inputLine, ';');
    for (const auto& pair : keyValuePairs) {
        std::vector<std::string> kv = splitString(pair, '=');
        if (kv.size() == 2) {
            jsonData[kv[0]] = processKeyValue1(kv[1]);
        } else {
          std::cerr << "Unexpected key=value data: " << pair << std::endl;
        }
    }
}

// Function to parse the ClinVar VCF file
std::unordered_map<std::string, ClinVarData> parseClinVarVCF(const std::string& filename, const std::regex& geneFilter) {
    std::unordered_map<std::string, ClinVarData> clinvarRecords;
    std::ifstream vcfFile(filename);
    std::string line;

    if (!vcfFile.is_open()) {
        std::cerr << "Error opening ClinVar VCF file!" << std::endl;
        return clinvarRecords;
    }

    while (getline(vcfFile, line)) {
        // Skip header lines
        if (line[0] == '#') continue;

        // Check if the line contains GENEINFO and match with the filter
        size_t geneInfoPos = line.find("GENEINFO=");
        if (geneInfoPos != std::string::npos) {
            size_t geneInfoEnd = line.find(';', geneInfoPos);
            if (geneInfoEnd == std::string::npos) { geneInfoEnd = line.length(); }

            std::string geneInfo = line.substr(geneInfoPos + 9, geneInfoEnd - geneInfoPos - 9);


            if (std::regex_match(geneInfo, geneFilter)) {
                // Extract RS ID
                size_t rsPos = line.find("RS=");
                if (rsPos != std::string::npos) {
                    size_t rsEnd = line.find(';', rsPos);
                    if (rsEnd == std::string::npos) { rsEnd = line.length(); }
                    std::string rsIDpart = line.substr(rsPos + 3, rsEnd - rsPos - 3);
                    // Split the RS numbers by ':'
                    std::vector<std::string> rsIDs = splitString(rsIDpart, ':');
                    // Extract REF and ALT alleles
                    std::istringstream iss(line);
                    std::string chrom, pos, id, ref, alt, qual, filter, info;
                    iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info;

                    // Store the data for each RS ID
                    for (const auto& rs : rsIDs) {
                        clinvarRecords[rs] = { chrom, pos, id, ref, alt, qual, filter, geneInfo, info};
                    }
                }
            }
        }
    }
    return clinvarRecords;
}

// Function to process the AncestryDNA file and match against ClinVar data
void processAncestryDNA(const std::string& filename, const std::unordered_map<std::string, ClinVarData>& clinvarRecords, bool onlyDiscrepancies, bool unmatchedRsId) {
    std::ifstream dnaFile(filename);
    std::string line;
    int processedLines;
    json jsonList;
    json jsonData;
    if (!dnaFile.is_open()) {
        std::cerr << "Error opening AncestryDNA file!" << std::endl;
        return;
    }
    processedLines = -1;
    while (getline(dnaFile, line)) {
        if (line[0] == '#') continue;  // Skip comment lines
        if (++processedLines==0) continue; // Skip first line

        trim_endlines(line);
        std::istringstream iss(line);
        std::string rsid, chromosome, position, allele1, allele2;
        iss >> rsid >> chromosome >> position >> allele1 >> allele2;
        if (rsid.substr(0, 2) == "rs") rsid = rsid.substr(2);

        // Find RS ID in ClinVar records
        auto it = clinvarRecords.find(rsid);
        if (it != clinvarRecords.end()) {
            const ClinVarData& clinvarData = it->second; // get value
            jsonData = json();

            // Check allele discrepancies
            std::string prefix = "-/-";  // Default to no discrepancies
            if (allele1 != clinvarData.refAllele || allele2 != clinvarData.refAllele) {
                if (allele1 != clinvarData.refAllele && allele2 != clinvarData.refAllele) {
                    prefix = "+/+";
                } else {
                    prefix = "+/-";
                }

                // Output discrepancies
                jsonData["rsID"]       = rsid;
                jsonData["allele1"]    = allele1;
                jsonData["allele2"]    = allele2;
                jsonData["chromosome"] = chromosome;
                jsonData["position"]   = position;
                jsonData["REF"]    = clinvarData.refAllele;
                jsonData["ALT"]    = clinvarData.altAllele;
                jsonData["FILTER"] = clinvarData.filter;
                jsonData["QUAL"] = clinvarData.qual;
                jsonData["ID"] = clinvarData.id;
                if (clinvarData.pos != position) {
                  std::cerr << "Warning: clinvar position does not match dna raw data position: " << clinvarData.pos << " " << position << std::endl;
                  // std::exit(-1);
                }
                processKeyValue0(jsonData, clinvarData.info);
                jsonList += jsonData;
                std::cout << prefix << "\t" << line << "\t" << clinvarData.info << std::endl;

            } else {
                if (!onlyDiscrepancies) {
                    // Output all matches if not filtering discrepancies
                    jsonData["rsID"]       = rsid;
                    jsonData["allele1"]    = allele1;
                    jsonData["allele2"]    = allele2;
                    jsonData["chromosome"] = chromosome;
                    jsonData["position"]   = position;
                    jsonData["REF"]    = clinvarData.refAllele;
                    jsonData["ALT"]    = clinvarData.altAllele;
                    jsonData["FILTER"] = clinvarData.filter;
                    jsonData["QUAL"] = clinvarData.qual;
                    jsonData["ID"] = clinvarData.id;
                    if (clinvarData.pos != position) {
                      std::cerr << "Warning: clinvar position does not match dna raw data position: " << clinvarData.pos << " " << position << std::endl;
                    }
                    processKeyValue0(jsonData, clinvarData.info);
                    jsonList += jsonData;
                    std::cout << prefix << "\t" << line << "\t" << clinvarData.info << std::endl;
                }
            }
        } else {
          if (unmatchedRsId) {
            std::cout << "?/?\t" << line << "\tNo match in clinvar.vcf!" <<  std::endl;
          }
        }
    }
    std::ofstream outfile("filtered-DNA.json");
    if (!outfile.is_open()) {
      std::cerr << "Error: Can't open filtered-DNA.json for write." << std::endl;
      std::exit(-1);
    } else {
      outfile << jsonList.dump();
    }
}

// Main function
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <clinvar.vcf> <AncestryDNA.txt> [options]" << std::endl;
        return 1;
    }

    std::string clinvarFile = argv[1];
    std::string ancestryFile = argv[2];

    // Process options (gene filtering and discrepancy output)
    std::regex geneFilter(".*");  // Default to no filter (match all genes)
    bool onlyDiscrepancies = false;
    bool unmatchedRsId = false;

    for (int i = 3; i < argc; ++i) {
        std::string option = argv[i];
        if (option == "--gene") {
            if (i + 1 < argc) {
                geneFilter = std::regex(argv[++i]);
            }
        } else if (option == "--discrepancies") {
            onlyDiscrepancies = true;
        } else if (option == "--unmatched-rsid") {
            unmatchedRsId = true;
        }
    }

    // Parse the ClinVar VCF file
    auto clinvarRecords = parseClinVarVCF(clinvarFile, geneFilter);

    // Process the AncestryDNA file and match with ClinVar data
    processAncestryDNA(ancestryFile, clinvarRecords, onlyDiscrepancies, unmatchedRsId);

    return 0;
}
