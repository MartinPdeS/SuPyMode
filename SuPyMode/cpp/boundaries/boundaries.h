#pragma once

#include <string>
#include <stdexcept>


enum class BoundaryCondition {
    Zero,
    Symmetric,
    AntiSymmetric
};


static constexpr std::string_view to_string(const BoundaryCondition bc) {
    switch (bc) {
        case BoundaryCondition::Zero:          return "Zero";
        case BoundaryCondition::Symmetric:     return "Symmetric";
        case BoundaryCondition::AntiSymmetric: return "AntiSymmetric";
    }
}

class Boundaries {
public:
    BoundaryCondition left;
    BoundaryCondition right;
    BoundaryCondition top;
    BoundaryCondition bottom;

    Boundaries() = default;

    Boundaries(
        const BoundaryCondition& left,
        const BoundaryCondition& right,
        const BoundaryCondition& top,
        const BoundaryCondition& bottom
    ) : left(left), right(right), top(top), bottom(bottom) {
        assert_valid();
    }

    Boundaries(
        const std::string& left,
        const std::string& right,
        const std::string& top,
        const std::string& bottom
    ) {
        this->left = parse_boundary(left);
        this->right = parse_boundary(right);
        this->top = parse_boundary(top);
        this->bottom = parse_boundary(bottom);
        assert_valid();
    }

    ~Boundaries() = default;

    BoundaryCondition parse_boundary(const std::string& bc_str) const {
        if (bc_str == "zero") return BoundaryCondition::Zero;
        if (bc_str == "symmetric") return BoundaryCondition::Symmetric;
        if (bc_str == "anti-symmetric") return BoundaryCondition::AntiSymmetric;
        throw std::invalid_argument("Boundary conditions must be 'zero', 'symmetric', or 'anti-symmetric'.");
    };


    void assert_valid() const {
        for (const BoundaryCondition& bc : {left, right, top, bottom}) {
            if (bc != BoundaryCondition::Zero && bc != BoundaryCondition::Symmetric && bc != BoundaryCondition::AntiSymmetric) {
                throw std::invalid_argument("Boundary conditions must be 'zero', 'symmetric', or 'anti-symmetric'.");
            }
        }
    }

    bool operator==(const Boundaries& other) const {
        return (
            left == other.left &&
            right == other.right &&
            top == other.top &&
            bottom == other.bottom
        );
    }

    void print() const {
        printf(
            "Boundaries:\n"
            "\tLeft: %.*s\n"
            "\tRight: %.*s\n"
            "\tTop: %.*s\n"
            "\tBottom: %.*s\n",
            static_cast<int>(to_string(left).size()),   to_string(left).data(),
            static_cast<int>(to_string(right).size()),  to_string(right).data(),
            static_cast<int>(to_string(top).size()),    to_string(top).data(),
            static_cast<int>(to_string(bottom).size()), to_string(bottom).data()
        );
    }




};
