function rootDir = generateRootDir(drctry)
    mkdir(['../plots/' datestr(today,'yyyy-mm-dd') '/' drctry]);
    rootDir=['../plots/' datestr(today,'yyyy-mm-dd') '/' drctry];
end