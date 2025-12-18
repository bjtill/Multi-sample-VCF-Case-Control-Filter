//BT Oct 24, 2025
//g++ -o casecontrol CaseControlVCFFilter_SingleInput_V2.cpp -std=c++11
//./vcf_filter_single --vcf input.vcf.gz --samples samples.txt --output filtered.vcf --min-depth 10 --max-depth 200 --percent-unique 80 --multi-allelic

//Single VCF filtering program using case-control filtering principles
//Filters variants based on: BED regions (optional), depth, genotype validity, and uniqueness
//Memory-efficient: Uses streaming decompression for .gz files to handle large (8GB+) VCF files

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <cstdio>

// Structure to hold BED region
struct BedRegion {
    std::string chrom;
    int start;
    int end;
    std::string gene;
};

// Structure to hold sample status
struct SampleStatus {
    char status; // '+', '-', or '0'
};

class VCFFilterSingle {
private:
    std::string vcfFile;
    std::string bedFile;
    std::string sampleFile;
    std::string outputFile;
    
    // Optional parameters
    int minDepth;
    int maxDepth;
    bool allowMultiAllelic;
    bool allowMissing;
    double percentUniqueness;
    
    // Data structures
    std::vector<BedRegion> bedRegions;
    std::map<std::string, SampleStatus> sampleStatuses;
    std::ofstream logFile;
    
    // Statistics
    int totalVariants;
    int passedVariants;
    
public:
    VCFFilterSingle() : minDepth(-1), maxDepth(-1), allowMultiAllelic(false), 
                         allowMissing(false), percentUniqueness(100.0), totalVariants(0), passedVariants(0) {}
    
    void setParameters(const std::string& vcf, const std::string& bed, 
                      const std::string& sample, const std::string& output) {
        vcfFile = vcf;
        bedFile = bed;
        sampleFile = sample;
        outputFile = output;
    }
    
    void setOptionalParams(int minD, int maxD, bool multiAllelic, bool allowMiss, double pctUniq) {
        minDepth = minD;
        maxDepth = maxD;
        allowMultiAllelic = multiAllelic;
        allowMissing = allowMiss;
        percentUniqueness = pctUniq;
    }
    
    bool run() {
        // Open log file
        std::string logFileName = outputFile + ".log";
        logFile.open(logFileName);
        if (!logFile.is_open()) {
            std::cerr << "Error: Cannot create log file" << std::endl;
            return false;
        }
        
        logTime("Starting VCF filtering pipeline");
        logTime("Input VCF: " + vcfFile);
        logTime("Output VCF: " + outputFile);
        
        // Log parameters
        logParameters();
        
        // Load BED file if provided
        if (!bedFile.empty()) {
            if (!loadBedFile()) {
                logTime("Error loading BED file");
                return false;
            }
            logTime("Loaded " + std::to_string(bedRegions.size()) + " regions from BED file");
        } else {
            logTime("No BED file provided - all regions will be processed");
        }
        
        // Load sample file
        if (!loadSampleFile()) {
            logTime("Error loading sample file");
            return false;
        }
        logTime("Loaded " + std::to_string(sampleStatuses.size()) + " sample statuses");
        
        // Count samples by status
        int numCases = 0, numControls = 0, numUnknown = 0;
        for (const auto& pair : sampleStatuses) {
            if (pair.second.status == '+') numCases++;
            else if (pair.second.status == '-') numControls++;
            else numUnknown++;
        }
        logTime("Cases: " + std::to_string(numCases) + 
                ", Controls: " + std::to_string(numControls) + 
                ", Unknown: " + std::to_string(numUnknown));
        
        // Process VCF file
        if (!processVCF()) {
            logTime("Error processing VCF file");
            return false;
        }
        
        logTime("Total variants processed: " + std::to_string(totalVariants));
        logTime("Variants passing filters: " + std::to_string(passedVariants));
        if (totalVariants > 0) {
            double pct = (double)passedVariants / totalVariants * 100.0;
            logTime("Retention rate: " + std::to_string(pct) + "%");
        }
        
        logTime("Pipeline completed successfully");
        logFile.close();
        
        return true;
    }
    
private:
    void logTime(const std::string& message) {
        time_t now = time(0);
        char buf[80];
        strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", localtime(&now));
        logFile << "[" << buf << "] " << message << std::endl;
        std::cout << "[" << buf << "] " << message << std::endl;
    }
    
    void logParameters() {
        logFile << "\n=== Filtering Parameters ===" << std::endl;
        if (!bedFile.empty()) {
            logFile << "BED file: " << bedFile << std::endl;
        } else {
            logFile << "BED file: None (all regions)" << std::endl;
        }
        logFile << "Sample file: " << sampleFile << std::endl;
        
        if (minDepth > 0) {
            logFile << "Minimum depth: " << minDepth << std::endl;
        } else {
            logFile << "Minimum depth: None" << std::endl;
        }
        
        if (maxDepth > 0) {
            logFile << "Maximum depth: " << maxDepth << std::endl;
        } else {
            logFile << "Maximum depth: None" << std::endl;
        }
        
        logFile << "Allow multi-allelic: " << (allowMultiAllelic ? "Yes" : "No") << std::endl;
        logFile << "Allow missing data: " << (allowMissing ? "Yes" : "No") << std::endl;
        logFile << "Percent uniqueness: " << percentUniqueness << "%" << std::endl;
        logFile << "===========================\n" << std::endl;
    }
    
    bool loadBedFile() {
        if (bedFile.empty()) {
            return true;
        }
        
        std::ifstream file(bedFile);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open BED file: " << bedFile << std::endl;
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            BedRegion region;
            iss >> region.chrom >> region.start >> region.end;
            
            // Check for optional 4th column (gene name)
            if (iss >> region.gene) {
                // Gene name present
            } else {
                region.gene = ".";
            }
            
            bedRegions.push_back(region);
        }
        
        file.close();
        return true;
    }
    
    bool loadSampleFile() {
        std::ifstream file(sampleFile);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open sample file: " << sampleFile << std::endl;
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            std::istringstream iss(line);
            std::string sample;
            char status;
            
            if (iss >> sample >> status) {
                SampleStatus ss;
                ss.status = status;
                sampleStatuses[sample] = ss;
            }
        }
        
        file.close();
        return true;
    }
    
    bool isInBedRegion(const std::string& chrom, int pos) {
        // If no BED regions loaded, accept all positions
        if (bedRegions.empty()) {
            return true;
        }
        
        for (const auto& region : bedRegions) {
            if (region.chrom == chrom && pos >= region.start && pos <= region.end) {
                return true;
            }
        }
        return false;
    }
    
    // Check if file is gzipped
    bool isGzipped(const std::string& filename) {
        return (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz");
    }
    
    // Open file with streaming decompression if needed
    FILE* openVCFFile(const std::string& filename) {
        if (isGzipped(filename)) {
            // Use popen to stream gunzip output - no temp file needed!
            std::string cmd = "gunzip -c " + filename;
            return popen(cmd.c_str(), "r");
        } else {
            return fopen(filename.c_str(), "r");
        }
    }
    
    void closeVCFFile(FILE* fp, bool isGzipped) {
        if (isGzipped) {
            pclose(fp);
        } else {
            fclose(fp);
        }
    }
    
    bool processVCF() {
        // Open file with streaming decompression if needed
        bool isGz = isGzipped(vcfFile);
        FILE* inFile = openVCFFile(vcfFile);
        
        if (!inFile) {
            logTime("Error: Cannot open VCF file: " + vcfFile);
            return false;
        }
        
        std::ofstream outFile(outputFile);
        if (!outFile.is_open()) {
            logTime("Error: Cannot create output file: " + outputFile);
            closeVCFFile(inFile, isGz);
            return false;
        }
        
        char* lineBuffer = new char[1048576]; // 1MB buffer for very long lines
        std::vector<std::string> sampleNames;
        int formatIdx = -1;
        bool headerProcessed = false;
        
        // Process file line by line using C-style file reading for better performance
        while (fgets(lineBuffer, 1048576, inFile) != NULL) {
            std::string line(lineBuffer);
            
            // Remove trailing newline
            if (!line.empty() && line[line.length()-1] == '\n') {
                line.erase(line.length()-1);
            }
            if (!line.empty() && line[line.length()-1] == '\r') {
                line.erase(line.length()-1);
            }
            
            // Pass through meta-information and header lines
            if (line[0] == '#') {
                outFile << line << std::endl;
                
                // Parse sample names from header line
                if (line.substr(0, 6) == "#CHROM") {
                    std::istringstream iss(line);
                    std::string token;
                    int idx = 0;
                    
                    while (iss >> token) {
                        if (idx >= 9) { // Sample names start at column 9
                            sampleNames.push_back(token);
                        }
                        if (token == "FORMAT") {
                            formatIdx = idx;
                        }
                        idx++;
                    }
                    headerProcessed = true;
                }
                continue;
            }
            
            if (!headerProcessed) {
                outFile << line << std::endl;
                continue;
            }
            
            totalVariants++;
            
            // Log progress every 10000 variants
            if (totalVariants % 10000 == 0) {
                logTime("Processed " + std::to_string(totalVariants) + " variants...");
            }
            
            // Parse variant line
            std::istringstream iss(line);
            std::string chrom, id, ref, alt, qual, filter, info, format;
            int pos;
            
            iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
            
            // Check BED region filter
            if (!isInBedRegion(chrom, pos)) {
                continue;
            }
            
            // Parse sample data
            std::vector<std::string> sampleData;
            std::string sample;
            while (iss >> sample) {
                sampleData.push_back(sample);
            }
            
            // Get format field indices
            std::vector<std::string> formatFields;
            std::istringstream formatIss(format);
            std::string field;
            while (std::getline(formatIss, field, ':')) {
                formatFields.push_back(field);
            }
            
            int gtIdx = -1, dpIdx = -1;
            for (size_t i = 0; i < formatFields.size(); i++) {
                if (formatFields[i] == "GT") gtIdx = i;
                if (formatFields[i] == "DP") dpIdx = i;
            }
            
            if (gtIdx == -1) {
                continue; // No GT field, skip
            }
            
            // Build genotype and depth maps
            std::map<std::string, std::string> genotypes;
            std::map<std::string, int> depths;
            
            for (size_t i = 0; i < sampleNames.size() && i < sampleData.size(); i++) {
                std::vector<std::string> fields;
                std::istringstream fieldIss(sampleData[i]);
                std::string f;
                while (std::getline(fieldIss, f, ':')) {
                    fields.push_back(f);
                }
                
                if (gtIdx < (int)fields.size()) {
                    genotypes[sampleNames[i]] = fields[gtIdx];
                }
                
                if (dpIdx != -1 && dpIdx < (int)fields.size()) {
                    try {
                        depths[sampleNames[i]] = std::stoi(fields[dpIdx]);
                    } catch (...) {
                        depths[sampleNames[i]] = 0;
                    }
                }
            }
            
            // Apply filters
            if (!passesDepthFilter(depths)) {
                continue;
            }
            
            if (!passesGenotypeValidityFilter(genotypes)) {
                continue;
            }
            
            // Check if at least one case and one control have non-missing genotypes
            if (!hasNonMissingCaseAndControl(genotypes)) {
                continue;
            }
            
            if (!passesUniquenessFilter(genotypes)) {
                continue;
            }
            
            // Variant passed all filters
            passedVariants++;
            outFile << line << std::endl;
        }
        
        delete[] lineBuffer;
        closeVCFFile(inFile, isGz);
        outFile.close();
        
        return true;
    }
    
    bool passesDepthFilter(const std::map<std::string, int>& depths) {
        if (minDepth <= 0 && maxDepth <= 0) {
            return true; // No depth filtering
        }
        
        for (const auto& pair : depths) {
            std::string sample = pair.first;
            int dp = pair.second;
            
            // Only check depth for samples with known status
            if (sampleStatuses.find(sample) == sampleStatuses.end()) {
                continue;
            }
            
            char status = sampleStatuses[sample].status;
            if (status == '+' || status == '-') {
                if (minDepth > 0 && dp < minDepth) return false;
                if (maxDepth > 0 && dp > maxDepth) return false;
            }
        }
        
        return true;
    }
    
    bool passesGenotypeValidityFilter(const std::map<std::string, std::string>& genotypes) {
        for (const auto& pair : genotypes) {
            std::string sample = pair.first;
            std::string gt = pair.second;
            
            // Only check genotypes for samples with known status
            if (sampleStatuses.find(sample) == sampleStatuses.end()) {
                continue;
            }
            
            char status = sampleStatuses[sample].status;
            if (status == '+' || status == '-') {
                if (!isValidGenotype(gt)) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    bool hasNonMissingCaseAndControl(const std::map<std::string, std::string>& genotypes) {
        bool hasNonMissingCase = false;
        bool hasNonMissingControl = false;
        
        for (const auto& pair : genotypes) {
            std::string sample = pair.first;
            std::string gt = pair.second;
            
            // Only check samples with known status
            if (sampleStatuses.find(sample) == sampleStatuses.end()) {
                continue;
            }
            
            char status = sampleStatuses[sample].status;
            
            // Check if genotype is non-missing
            bool isMissing = (gt == "./." || gt == ".|." || gt.empty() || gt == ".");
            
            if (!isMissing && isValidGenotype(gt)) {
                if (status == '+') {
                    hasNonMissingCase = true;
                } else if (status == '-') {
                    hasNonMissingControl = true;
                }
            }
            
            // Early exit if both found
            if (hasNonMissingCase && hasNonMissingControl) {
                return true;
            }
        }
        
        return hasNonMissingCase && hasNonMissingControl;
    }
    
    bool isValidGenotype(const std::string& gt) {
        if (gt.empty() || gt == ".") {
            return false;
        }
        
        // If allowing missing data, ./. is considered valid (neutral)
        if (allowMissing && (gt == "./." || gt == ".|.")) {
            return true;
        }
        
        // If not allowing missing data, ./. is invalid
        if (!allowMissing && (gt == "./." || gt == ".|.")) {
            return false;
        }
        
        if (allowMultiAllelic) {
            // Allow any numeric genotype
            return true;
        } else {
            // Only allow 0/0, 0/1, 1/0, 1/1
            return (gt == "0/0" || gt == "0/1" || gt == "1/0" || gt == "1/1" ||
                    gt == "0|0" || gt == "0|1" || gt == "1|0" || gt == "1|1");
        }
    }
    
    std::set<std::string> extractAlleles(const std::string& gt) {
        std::set<std::string> alleles;
        std::string allele;
        for (char c : gt) {
            if (c == '/' || c == '|') {
                if (!allele.empty()) alleles.insert(allele);
                allele.clear();
            } else if (c != '.') {
                allele += c;
            }
        }
        if (!allele.empty()) alleles.insert(allele);
        return alleles;
    }
    
    bool passesUniquenessFilter(const std::map<std::string, std::string>& genotypes) {
        // Collect all case genotypes (skip missing if needed)
        std::set<std::string> caseGenotypes;
        int totalCases = 0;
        
        for (const auto& pair : genotypes) {
            std::string sample = pair.first;
            std::string gt = pair.second;
            
            char status = '0';
            if (sampleStatuses.find(sample) != sampleStatuses.end()) {
                status = sampleStatuses[sample].status;
            }
            
            if (status == '+') {
                if (!isValidGenotype(gt)) continue;
                
                // Skip missing genotypes when allowing missing data
                if (allowMissing && (gt == "./." || gt == ".|.")) continue;
                
                // If not allowing missing, count it but skip collecting genotype
                if (!allowMissing && (gt == "./." || gt == ".|.")) {
                    totalCases++;
                    continue;
                }
                
                totalCases++;
                std::string normalizedGt = normalizeGenotype(gt);
                caseGenotypes.insert(normalizedGt);
            }
        }
        
        // Collect all control genotypes (skip missing if needed)
        std::set<std::string> controlGenotypes;
        int totalControls = 0;
        
        for (const auto& pair : genotypes) {
            std::string sample = pair.first;
            std::string gt = pair.second;
            
            char status = '0';
            if (sampleStatuses.find(sample) != sampleStatuses.end()) {
                status = sampleStatuses[sample].status;
            }
            
            if (status == '-') {
                if (!isValidGenotype(gt)) continue;
                
                // Skip missing genotypes when allowing missing data
                if (allowMissing && (gt == "./." || gt == ".|.")) continue;
                
                // If not allowing missing, count it but skip collecting genotype
                if (!allowMissing && (gt == "./." || gt == ".|.")) {
                    totalControls++;
                    continue;
                }
                
                totalControls++;
                std::string normalizedGt = normalizeGenotype(gt);
                controlGenotypes.insert(normalizedGt);
            }
        }
        
        if (totalCases == 0 || totalControls == 0) return false;
        
        // Check cases for uniqueness (different from controls)
        int uniqueCases = 0;
        
        for (const auto& pair : genotypes) {
            std::string sample = pair.first;
            std::string gt = pair.second;
            
            char status = '0';
            if (sampleStatuses.find(sample) != sampleStatuses.end()) {
                status = sampleStatuses[sample].status;
            }
            
            if (status == '+') {
                if (!isValidGenotype(gt)) continue;
                
                // If allowMissing is true, skip missing genotypes
                if (allowMissing && (gt == "./." || gt == ".|.")) continue;
                
                // If not allowing missing, missing genotype is considered unique
                if (!allowMissing && (gt == "./." || gt == ".|.")) {
                    uniqueCases++;
                    continue;
                }
                
                // Normalize the case genotype
                std::string normalizedGt = normalizeGenotype(gt);
                
                // Check if this case genotype is DIFFERENT from all control genotypes
                if (controlGenotypes.find(normalizedGt) == controlGenotypes.end()) {
                    uniqueCases++;
                }
            }
        }
        
        // Check controls for uniqueness (different from cases)
        int uniqueControls = 0;
        
        for (const auto& pair : genotypes) {
            std::string sample = pair.first;
            std::string gt = pair.second;
            
            char status = '0';
            if (sampleStatuses.find(sample) != sampleStatuses.end()) {
                status = sampleStatuses[sample].status;
            }
            
            if (status == '-') {
                if (!isValidGenotype(gt)) continue;
                
                // If allowMissing is true, skip missing genotypes
                if (allowMissing && (gt == "./." || gt == ".|.")) continue;
                
                // If not allowing missing, missing genotype is considered unique
                if (!allowMissing && (gt == "./." || gt == ".|.")) {
                    uniqueControls++;
                    continue;
                }
                
                // Normalize the control genotype
                std::string normalizedGt = normalizeGenotype(gt);
                
                // Check if this control genotype is DIFFERENT from all case genotypes
                if (caseGenotypes.find(normalizedGt) == caseGenotypes.end()) {
                    uniqueControls++;
                }
            }
        }
        
        // Calculate percentages
        double casePct = (double)uniqueCases / totalCases * 100.0;
        double controlPct = (double)uniqueControls / totalControls * 100.0;
        
        // BOTH cases and controls must meet the threshold
        return (casePct >= percentUniqueness && controlPct >= percentUniqueness);
    }
    
    // Helper function to normalize genotypes (e.g., 0/1 -> 0/1, 1/0 -> 0/1, 0|1 -> 0/1)
    // Phased (|) and unphased (/) are treated as identical
    std::string normalizeGenotype(const std::string& gt) {
        // Extract alleles
        std::vector<std::string> alleles;
        std::string allele;
        
        for (char c : gt) {
            if (c == '/' || c == '|') {
                if (!allele.empty()) {
                    alleles.push_back(allele);
                    allele.clear();
                }
            } else if (c != '.') {
                allele += c;
            }
        }
        if (!allele.empty()) {
            alleles.push_back(allele);
        }
        
        // Sort alleles to normalize (0/1 and 1/0 become the same)
        if (alleles.size() == 2) {
            std::sort(alleles.begin(), alleles.end());
            // Always use '/' separator (ignore phasing)
            return alleles[0] + "/" + alleles[1];
        } else if (alleles.size() == 1) {
            return alleles[0] + "/" + alleles[0];
        }
        
        return gt;
    }
};

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " [OPTIONS]" << std::endl;
    std::cout << "\nRequired arguments:" << std::endl;
    std::cout << "  --vcf <file>           Input multi-sample VCF file (can be .vcf or .vcf.gz)" << std::endl;
    std::cout << "  --samples <file>       Sample status file (sample_name status per line)" << std::endl;
    std::cout << "  -o, --output <file>    Output VCF file" << std::endl;
    std::cout << "\nOptional arguments:" << std::endl;
    std::cout << "  --bed <file>           BED file with regions (optional, no BED = all regions)" << std::endl;
    std::cout << "  --min-depth <int>      Minimum depth (DP) for case/control samples" << std::endl;
    std::cout << "  --max-depth <int>      Maximum depth (DP) for case/control samples" << std::endl;
    std::cout << "  --multi-allelic        Allow multi-allelic variants (0/2, 1/2, etc.)" << std::endl;
    std::cout << "  -a, --allow-missing    Allow missing data (./. treated as neutral)" << std::endl;
    std::cout << "  --percent-unique <pct> Percent uniqueness threshold (default: 100)" << std::endl;
    std::cout << "  -h, --help             Show this help message" << std::endl;
    std::cout << "\nSample status file format:" << std::endl;
    std::cout << "  sample1 +   (case)" << std::endl;
    std::cout << "  sample2 -   (control)" << std::endl;
    std::cout << "  sample3 0   (unknown/ignore)" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }
    
    std::string vcfFile, bedFile, sampleFile, outputFile;
    int minDepth = -1, maxDepth = -1;
    bool multiAllelic = false;
    bool allowMissing = false;
    double percentUnique = 100.0;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--vcf" && i + 1 < argc) {
            vcfFile = argv[++i];
        } else if (arg == "--bed" && i + 1 < argc) {
            bedFile = argv[++i];
        } else if (arg == "--samples" && i + 1 < argc) {
            sampleFile = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (arg == "--min-depth" && i + 1 < argc) {
            minDepth = std::stoi(argv[++i]);
        } else if (arg == "--max-depth" && i + 1 < argc) {
            maxDepth = std::stoi(argv[++i]);
        } else if (arg == "--multi-allelic") {
            multiAllelic = true;
        } else if (arg == "-a" || arg == "--allow-missing") {
            allowMissing = true;
        } else if (arg == "--percent-unique" && i + 1 < argc) {
            percentUnique = std::stod(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        }
    }
    
    // Validate required arguments
    if (vcfFile.empty() || sampleFile.empty() || outputFile.empty()) {
        std::cerr << "Error: Missing required arguments" << std::endl;
        printUsage(argv[0]);
        return 1;
    }
    
    // Run the filter
    VCFFilterSingle filter;
    filter.setParameters(vcfFile, bedFile, sampleFile, outputFile);
    filter.setOptionalParams(minDepth, maxDepth, multiAllelic, allowMissing, percentUnique);
    
    if (!filter.run()) {
        std::cerr << "Error: Pipeline failed" << std::endl;
        return 1;
    }
    
    std::cout << "\nPipeline completed successfully!" << std::endl;
    std::cout << "Output VCF: " << outputFile << std::endl;
    std::cout << "Log file: " << outputFile << ".log" << std::endl;
    
    return 0;
}
