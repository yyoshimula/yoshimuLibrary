function results = runAllTests()
% runAllTests  Run the yoshimuLibrary test suite.
%
%   results = runAllTests() adds the library to the path, runs every test
%   file directly inside the tests/ folder (NON-recursively, so the
%   interactive verification scripts in tests/manual are never executed),
%   displays a summary table, and returns the TestResult array.
%
%   Note: runtests(folder) is non-recursive by default ('IncludeSubfolders'
%   defaults to false), which keeps tests/manual out of the suite. Test
%   files are still selected explicitly below so that this function itself
%   (whose name ends in "Tests") is never mistaken for a test file.

root = fileparts(fileparts(mfilename('fullpath')));
% genpath includes hidden folders such as .claude/worktrees, whose stale
% library copies would shadow the real files
p = strsplit(genpath(root), pathsep);
p = p(~cellfun(@isempty, p) & ~contains(p, [filesep '.']));
addpath(strjoin(p, pathsep));

testsDir = fullfile(root, 'tests');

% Select test files explicitly: *.m files directly in tests/ whose names
% start or end with "test" (MATLAB's test-file naming convention),
% excluding this runner itself.
listing = dir(fullfile(testsDir, '*.m'));
names = {listing.name};
keep = false(size(names));
for i = 1:numel(names)
    [~, base] = fileparts(names{i});
    isTestName = ~isempty(regexpi(base, '^test', 'once')) || ...
        ~isempty(regexpi(base, 'test[s]?$', 'once'));
    keep(i) = isTestName && ~strcmpi(base, 'runAllTests');
end
testFiles = fullfile(testsDir, names(keep));

if isempty(testFiles)
    error('runAllTests:noTests', 'No test files found in %s.', testsDir);
end

results = runtests(testFiles);

disp(table(results));

end
