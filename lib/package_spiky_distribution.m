function package_spiky_distribution
% PACKAGE_SPIKY_DISTRIBUTION
% Packages the Spiky software in a zip file for distribution.
%
%  - depends on Linux environment to run
%  - zip file 
%  - include build number in package name
%  - exclude .asv, .mbk and .mat files
%  - exclude .svn directories
%  - exclude ./bin/src directory
%

sArchName = ['spiky_build_' spiky('GetBuildNumber') '.zip'];

sPath = which('spiky');
sPath = sPath(1:end-13);
cd(sPath)

% initial compress (all files)
sCmd = ['zip -r ' sArchName ' ./spiky/'];
system(sCmd)

% delete .asv files
sCmd = ['zip -d ' sArchName ' \*.asv']
system(sCmd)

% delete .mat files
sCmd = ['zip -d ' sArchName ' \*.mat']
system(sCmd)

% delete backup directories
sCmd = ['zip -d ' sArchName ' \*backup\*']
system(sCmd)

% delete .mbk files
sCmd = ['zip -d ' sArchName ' \*.mbk']
system(sCmd)

% delete .bak files
sCmd = ['zip -d ' sArchName ' \*.bak']
system(sCmd)

% delete .bak files
sCmd = ['zip -d ' sArchName ' \*.exe']
system(sCmd)

% delete .pdf files
sCmd = ['zip -d ' sArchName ' \*.pdf']
system(sCmd)

% delete .svn directories
sCmd = ['zip -d ' sArchName ' \*.svn\*']
system(sCmd)

% delete wiki directory
sCmd = ['zip -d ' sArchName ' \*wiki\*']
system(sCmd)

% delete .m~ files
sCmd = ['zip -d ' sArchName ' \*.m~']
system(sCmd)


