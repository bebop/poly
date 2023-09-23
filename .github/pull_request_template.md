## Changes in this PR
*Clearly and concisely summarize the changes you are making. Bullet points are completely okay. Please be specific, saying "improves X" is not enough!*

### Why are you making these changes?
*Explain why these changes are necessary. Link to GitHub issues here with the format `fixes: #XXX` to indicate this PR resolves the issue.*

### Are any changes breaking? (IMPORTANT)
*Will merging this PR change `poly`'s API in a non-backwards-compatible manner?*

*Examples of breaking changes:*
* *Removing a method from a struct.*
* *Deleting/moving a package.*
* *Adding a method to an interface (client may have made their own     implementation of this interface, and adding a method to the     interface could cause client's implementation to no longer satisfy     the interface).*

*Examples of non-breaking changes:*
* *Adding a method to a struct.*
* *Adding a function.*
* *Fixing a bug in a function/method.*
* *Creating a new package.*

## Pre-merge checklist
*All of these must be satisfied before this PR is considered
ready for merging. Mergeable PRs will be prioritized for review.*

* [ ] New packages/exported functions have docstrings.
* [ ] New/changed functionality is thoroughly tested.
* [ ] New/changed functionality has a function giving an example of its usage in the associated test file. See `primers/primers_test.go` for what this might look like.
* [ ] Changes are documented in `CHANGELOG.md` in the `[Unreleased]` section.
* [ ] All code is properly formatted and linted.
* [ ] The PR template is filled out.
