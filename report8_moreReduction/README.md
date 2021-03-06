# Report 7 *Further reduction*

1. [x] constructing the before75_after75 cDBG.
2. [ ] Construct a graph from the unitigs.
3. [ ] Extract all the free-end nodes from #2 to a new fasta file.
    - Fasta file headers should have the component serial ID including this node & Node ID.
    - Preserving the previous information will allow us to prevent the circular paths in free-ends within the same component. The extension should only be allowed with nodes from other component.
4. [ ] Construct a cDBG from the fasta file in #3 with k=25.
    - [ ] 1st scenario (most common): unconnected unitigs with len=75nt (edge that did not attach with any other node) SHOULD BE EXCLUDED.
    - [ ] 2nd scenario: Found that there are duplicate 25-mers within the 75nt unitigs which will cause multiple unitigs generation (Total len=75nt of the formed isolated component) SHOULD BE EXCLUDED [Circularization].
    - [ ] 3rd scenario: component length > 75nt and when searching in the Fasta file, we will find that **only two unitigs** aligned to a component, these two edges(unitigs) from two separate components should be merged to single component.
    - [ ] 4th scenario: If matched multiple unitigs in a single components, create a cDBG from the matched unitigs free end edges in the components by k=k+2=27.
5. [ ] Repeat #4.