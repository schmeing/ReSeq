from pathlib import Path

import DataStats


def main() -> None:
    print(Path(DataStats.__file__).name)


if __name__ == "__main__":
    main()
