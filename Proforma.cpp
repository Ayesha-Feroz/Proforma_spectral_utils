#include <regex>

std::regex MOD_GLOBAL_L("<");
std::regex MOD_GLOBAL_R(">");
std::regex MOD_L("\\[");
std::regex MOD_R("\\]");
std::regex MOD_LABILE_L("\\{");
std::regex MOD_LABILE_R("\\}");
std::regex MOD_RANGE_L("\\(");
std::regex MOD_RANGE_R("\\)");
std::regex AA("^[A-Z]$");
std::regex LETTER("^[A-Za-z]$");
std::regex INT("^\\d+$");
std::regex SIGNED_INT("^[+-]?\\d+$");
std::regex NUMBER("^(\\d+)?\\.?\\d*$");
std::regex SIGNED_NUMBER("^[+-]?(\\d+)?\\.?\\d*$");
std::regex WS("\\s+");
std::regex TEXT(".+");
std::regex MONOSACCHARIDE("[A-Z][a-z]?(?:-?[A-Za-z]+)*");
std::regex MOD_COUNT("\\d+");
std::regex CV_ABBREV_OPT("[UM]", std::regex_constants::icase);
std::regex CV_ABBREV("[RGX]", std::regex_constants::icase);
std::regex CV_NAME("(UNIMOD|MOD|RESID|XLMOD|GNO)", std::regex_constants::icase);
std::regex MOD_MASS_OBS("Obs", std::regex_constants::icase);
std::regex MOD_MASS("[+-]?\\d*(\\.\\d+)?");
std::regex FORMULA("[A-Za-z]+[+-]?\\d*\\s*");
std::regex MOD_LABEL_XL("XL");
std::regex MOD_LABEL_BRANCH("BRANCH");
std::regex MOD_LABEL("[A-Za-z0-9]+");
std::regex MOD_SCORE("[+-]?\\d*(\\.\\d+)?");

std::regex ionRegex("(.*?,)?(.*)");
std::smatch match;

struct Monosaccharide {
    std::string name;
    int count;
};

struct Mod {
    std::string name;
    std::string accession;
    std::string mass;
    std::string formula;
    std::vector<Monosaccharide> glycans;
    std::string info;
    std::string label;
    std::string score;
};

struct AminoAcid {
    std::string letter;
    std::vector<Mod> mods;
};

struct ModRangePos {
    std::vector<AminoAcid> aminoAcids;
    std::vector<Mod> mods;
};

struct Peptide {
    std::vector<Mod> modsGlobal;
    std::vector<Mod> modsUnknownPos;
    std::vector<Mod> modsLabile;
    Mod modNTerm;
    std::vector<AminoAcid> aminoAcids;
    Mod modCTerm;
};

struct Proteoform {
    std::vector<Peptide> peptides;
    std::string charge;
    std::string ion1;
    std::string ion2;
};

void parseIon(std::string ionString, Proteoform& proteoform) {
    if (std::regex_search(ionString, match, ionRegex)) {
        if (match[1].str().length() > 0) {
            proteoform.ion1 = match[1].str().substr(0, match[1].str().length()-1);
        }
        proteoform.ion2 = match[2].str();
    } else {
        proteoform.ion2 = ionString;
    }
}

void parseGlycan(std::string glycanString, Mod& mod) {
    std::regex_token_iterator
