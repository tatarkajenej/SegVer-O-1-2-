#include "segver.cpp"

//////////
// MAIN //
//////////

int main(int argc, char* args[]) {

    std::vector<std::vector<CoordConfigFamily>> picture = {};

    picture.push_back({CoordConfigFamily(1,2, "m", "n", 0, {"m+1"}, {"n+1"}, {"(-m^2+m+4+mn+n) / 2"})});

    picture.push_back({picture[0][0].inductant("m-2", "n-3m^2+6m-2")});
    picture[1].push_back(picture[1][0].inductant("m", "n-2"));
    picture[1].push_back(picture[1][1].inductant("m", "n-1"));

    picture.push_back({picture[1][0].inductant("m-2", "n-12m+22")});
    picture[2].push_back(picture[2][0].inductant("m", "n-2"));

    picture.push_back({picture[2][0].inductant("m-2", "n-12m+22")});
    picture[3].push_back(picture[3][0].inductant("m", "n-2"));

    picture.push_back({picture[3][0].inductant("m-2", "n-12m+22")});
    picture[4].push_back(picture[4][0].inductant("m", "n-2"));

    picture.push_back({picture[4][0].inductant("m-2", "n-12m+22")});

    picture.push_back({picture[5][0].inductant("m-2", "n-12m+22")});

    #include "select.cpp"
}
