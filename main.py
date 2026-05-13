from treeBuilder import build_trees, visualize


def get_positive_integer(prompt):
    """
    Repeatedly asks the user for a positive integer.
    """
    while True:
        user_input = input(prompt).strip()

        try:
            value = int(user_input)

            if value <= 0:
                print("Please enter a number greater than 0.")
            else:
                return value

        except ValueError:
            print("Invalid input. Please enter a whole number.")


def get_preference():
    """
    Asks the user whether they prefer speed or accuracy.
    """
    while True:
        preference = input("Do you prefer speed or accuracy? ").strip().lower()

        if preference in ["speed", "fast", "time"]:
            return "speed"

        if preference in ["accuracy", "accurate"]:
            return "accuracy"

        print("Please type either 'speed' or 'accuracy'.")


def get_yes_or_no(prompt):
    """
    Repeatedly asks the user for yes or no.
    """
    while True:
        answer = input(prompt).strip().lower()

        if answer in ["yes", "y"]:
            return True

        if answer in ["no", "n"]:
            return False

        print("Please type yes or no.")


def classify_dataset(taxa_count, sequence_length):
    """
    Classifies a dataset based on taxa count and sequence length.
    """
    if taxa_count <= 50:
        taxa_size = "small"
    elif taxa_count <= 500:
        taxa_size = "medium"
    else:
        taxa_size = "large"

    if sequence_length <= 500:
        sequence_size = "short"
    elif sequence_length <= 5000:
        sequence_size = "medium"
    else:
        sequence_size = "long"

    return taxa_size, sequence_size


def recommend_algorithm(taxa_count, sequence_length, preference):
    """
    Recommends a phylogenetic tree-building algorithm.

    The recommendation rules are based on common behavior of the algorithms:
    - UPGMA is fast, but assumes equal evolutionary rates.
    - NJ is usually a strong speed/accuracy balance.
    - ML is usually accuracy-focused but slower.
    - MP can be useful but may become computationally expensive.
    """
    taxa_size, sequence_size = classify_dataset(taxa_count, sequence_length)

    if preference == "speed":
        if taxa_size == "small":
            algorithm = "UPGMA"
            reason = (
                "UPGMA is recommended because the dataset is small and the user "
                "prioritized speed. UPGMA is a fast distance-based method and is "
                "useful as a quick baseline."
            )
        else:
            algorithm = "NJ"
            reason = (
                "Neighbor Joining is recommended because it is generally fast while "
                "handling medium or large datasets better than UPGMA when evolutionary "
                "rates differ."
            )

    else:
        if taxa_size == "small":
            algorithm = "ML"
            reason = (
                "Maximum Likelihood is recommended because the dataset is small enough "
                "for a more computationally intensive method, and the user prioritized "
                "accuracy."
            )
        elif taxa_size == "medium":
            algorithm = "ML"
            reason = (
                "Maximum Likelihood is recommended for accuracy, although it may take "
                "longer than distance-based methods. If runtime becomes a problem, NJ "
                "is the next best practical option."
            )
        else:
            algorithm = "NJ"
            reason = (
                "Neighbor Joining is recommended because the dataset is large. ML may "
                "be more accuracy-focused, but it can become computationally expensive "
                "as the number of taxa increases."
            )

    return {
        "algorithm": algorithm,
        "reason": reason,
        "taxa_size": taxa_size,
        "sequence_size": sequence_size
    }


def get_ml_model():
    """
    Lets the user choose a Maximum Likelihood substitution model.

    These are common/simple model choices. The full IQ-TREE model list is larger,
    but this keeps the UI understandable.
    """
    models = {
        "1": "JC",
        "2": "HKY",
        "3": "GTR",
        "4": "GTR+G"
    }

    print("\nMaximum Likelihood model options:")
    print("1 - JC")
    print("2 - HKY")
    print("3 - GTR")
    print("4 - GTR+G")

    while True:
        choice = input("Choose an ML model number: ").strip()

        if choice in models:
            return models[choice]

        print("Invalid choice. Please choose 1, 2, 3, or 4.")


def print_recommendation(taxa_count, sequence_length, preference, recommendation):
    """
    Prints the recommendation summary.
    """
    print("\n==============================")
    print("PhyloBench Recommendation")
    print("==============================")
    print(f"Number of taxa/sequences: {taxa_count}")
    print(f"Sequence length: {sequence_length}")
    print(f"Taxa category: {recommendation['taxa_size']}")
    print(f"Sequence length category: {recommendation['sequence_size']}")
    print(f"User preference: {preference}")
    print("\nRecommended algorithm:")
    print(recommendation["algorithm"])
    print("\nReason:")
    print(recommendation["reason"])
    print("==============================")


def main():
    print("==============================")
    print("PhyloBench Text UI")
    print("==============================")
    print("This tool recommends a phylogenetic tree-building algorithm.")
    print("It can also build a tree from a Clustal-format sequence alignment file.\n")

    sequence_path = input(
        "Enter the path to the sequence alignment file in Clustal format: "
    ).strip()

    taxa_count = get_positive_integer("Enter the number of taxa/sequences: ")
    sequence_length = get_positive_integer("Enter the sequence length: ")
    preference = get_preference()

    recommendation = recommend_algorithm(
        taxa_count,
        sequence_length,
        preference
    )

    print_recommendation(
        taxa_count,
        sequence_length,
        preference,
        recommendation
    )

    build_tree_choice = get_yes_or_no(
        "\nWould you like to build the recommended tree now? "
    )

    if not build_tree_choice:
        print("\nTree building skipped.")
        return

    recommended_algorithm = recommendation["algorithm"]

    ml_model = "JC"
    ml_random_seed = None

    if recommended_algorithm == "ML":
        choose_model = get_yes_or_no(
            "Would you like to choose the ML substitution model? "
        )

        if choose_model:
            ml_model = get_ml_model()

        use_seed = get_yes_or_no(
            "Would you like to set a random seed for ML reproducibility? "
        )

        if use_seed:
            ml_random_seed = get_positive_integer("Enter random seed: ")

    print(f"\nBuilding tree using {recommended_algorithm}...")

    try:
        trees = build_trees(
            sequence_path,
            model="identity",
            ml_model=ml_model,
            ml_random_seed=ml_random_seed,
            desired_trees=[recommended_algorithm]
        )

        print(f"\nSuccessfully built {recommended_algorithm} tree.")

        show_tree = get_yes_or_no(
            "Would you like to visualize the tree? "
        )

        if show_tree:
            visualize(trees, [recommended_algorithm])

    except FileNotFoundError:
        print("\nError: The file path was not found.")
        print("Please check the sequence path and try again.")

    except Exception as error:
        print("\nAn error occurred while building the tree.")
        print("Possible causes:")
        print("- The sequence file is not in Clustal format.")
        print("- The alignment contains invalid sequence characters.")
        print("- The chosen model is not compatible with the data.")
        print(f"\nError details: {error}")


if __name__ == "__main__":
    main()