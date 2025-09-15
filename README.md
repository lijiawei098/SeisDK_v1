# SeisDK_v1
A data-driven framework to systematically test the dragon-king earthquake hypothesis

This MATLAB program is designed to systematically test the dragon-king earthquake hypothesis.

## 🧑‍💻 Authors

- **Jiawei Li**  (lijw3@sustech.edu.cn | lijw@pku.edu.cn | lijw@cea-igp.ac.cn)

First version completed in **June 2025**.

### The framework we developed for detecting dragon-king earthquakes comprises four main steps, each carefully designed to address this critical challenge:
- (1). Defining the Study Region:<br>
  The first step focuses on identifying clusters of seismicity using kernel density estimation (KDE). This technique ensures that the study domain reflects natural groupings of earthquakes rather than arbitrary boundaries. For further analysis, this step is pivotal for identifying a set of events that are likely to be associated with the same seismogenesis. Studying the Gutenberg-Richter distribution within coherent seismicity regions, rather than in aggregate, is crucial to uncover localized geophysical variations and distinct earthquake types like dragon-kings, which broad-scale analysis obscures.
- (2). Identifying Dragon-King Candidates:<br>
  Once the study regions are defined, potential dragon-king earthquake candidates are identified by analyzing the statistical distribution of seismicity in each region separately. This involves examining the tail behavior of the FMD in each region and searching for anomalies that deviate markedly from the expected patterns of the generating process.
- (3). Rigorous Statistical Testing:<br>
  The next step applies rigorous statistical methods to individually test whether the identified candidates truly represent dragon-king earthquakes. These tests are specifically designed to ensure that the observed anomalies are not due to random variations or errors in data processing but are statistically significant deviations indicative of unique underlying processes.
- (4). Validation of Dragon-King Events:<br>
  If dragon-king earthquakes are detected, the final step involves validating these events through further statistical analysis using block tests. This may also include analyzing the geological and geophysical context and physical mechanisms associated with the events, and comparing them to similar phenomena in other regions or contexts.

## 📖 References

If you use this code, please cite or refer to the following studies:

Li, J.#, Sornette, D.# (2025). Haicheng and Tangshan earthquakes as potential dragon-kings. Science China Earth Sciences (under review) or https://arxiv.org/abs/2504.21310
