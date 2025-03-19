#include "segver.cpp"

//////////
// MAIN //
//////////

int main(int argc, char* args[]) {

    std::vector<std::vector<CoordConfigFamily>> picture = {};

    picture.push_back({CoordConfigFamily(1,2, "m", "n", 0, {"m+1"}, {"n+1"}, {"(-m^2+m+3+mn+n) / 2"})});  // ugly
    picture[0].push_back(CoordConfigFamily(1,2, "m", "n", 0, {"m+1"}, {"n+1"}, {"(-m^2+m+4+mn+n) / 2"}));  // nice

    picture.push_back({picture[0][0].inductant("m", "(n-m-1)/2", picture[0][1])}); // makes sense if n = m+1 (4)
    picture[1].push_back(picture[0][0].inductant("m", "(n-m-3)/2", picture[0][1]));    // makes sense if n = m-1 (4)

    picture[1].push_back(picture[1][0].inductant("m", "n-4"));
    picture[1].push_back(picture[1][1].inductant("m", "n-4"));
    picture[1].push_back(picture[1][2].inductant("m", "n-4"));
    picture[1].push_back(picture[1][3].inductant("m", "n-4"));

    picture.push_back({picture[1][0].inductant("m-2", "n-6")});
    picture[2].push_back(picture[1][1].inductant("m-2", "n-6"));
    picture[2].push_back(picture[2][0].inductant("m", "n-4"));
    picture[2].push_back(picture[2][1].inductant("m", "n-4"));
    picture[2].push_back(picture[2][2].inductant("m", "n-4"));
    picture[2].push_back(picture[2][3].inductant("m", "n-4"));

    picture.push_back({picture[2][0].inductant("m-2", "n-6")});
    picture[3].push_back(picture[2][1].inductant("m-2", "n-6"));
    picture[3].push_back(picture[3][0].inductant("m", "n-4"));
    picture[3].push_back(picture[3][1].inductant("m", "n-4"));

    picture.push_back({picture[3][0].inductant("m-2", "n-6")});
    picture[4].push_back(picture[3][1].inductant("m-2", "n-6"));

    #include "select.cpp"
}
