int_t which1 = atoi(args[1]);
int_t which2 = atoi(args[2]);
int_t u = atoi(args[3]);
int_t v = atoi(args[4]);
std::string filename;
if (args[5]) {
    filename = args[5];
}
else {
    std::string problemName = args[0];
    int_t lastSlash = problemName.rfind("/");
    if (lastSlash == std::string::npos) {
        problemName = "";
    }
    else {
        problemName = problemName.substr(lastSlash+1) + "-";
    }
    filename = "cert/" + problemName + std::to_string(which1) + "-" + std::to_string(which2) + "-" + std::to_string(u) + "-" + std::to_string(v) + ".run";
    std::cout << "Using implicitly deduced certificate filename: " << filename << std::endl;
}

std::cout << picture[which1][which2] << std::endl;

time_point t1 = std::chrono::system_clock::now();
CoordConfig cc = picture[which1][which2].eval(u,v);
time_point t2 = std::chrono::system_clock::now();
std::cout << "The CoordConfig object built in " << timeElapsed(t1,t2) << "." << std::endl;
std::cout << cc << std::endl;
cc.getMonomials(true);
cc.getVariables(true);
std::cout << "matrix dimensions: " << cc.monomialCount << " x " << cc.columnCount << std::endl;
bool result = cc.computeWithCertificate(filename);
if (result){
    std::cout << "NON-DEFECTIVE" << std::endl;
}
else {
    std::cout << "DEFECTIVE" << std::endl;
}

return 0;
