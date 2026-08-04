#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace Spatial
{
    struct Bounds2
    {
        float minX = 0.0f;
        float minZ = 0.0f;
        float maxX = 0.0f;
        float maxZ = 0.0f;
    };

    inline bool Intersects(const Bounds2& a, const Bounds2& b)
    {
        return a.minX <= b.maxX && a.maxX >= b.minX &&
            a.minZ <= b.maxZ && a.maxZ >= b.minZ;
    }

    inline bool Contains(const Bounds2& outer, const Bounds2& inner)
    {
        return inner.minX >= outer.minX && inner.maxX <= outer.maxX &&
            inner.minZ >= outer.minZ && inner.maxZ <= outer.maxZ;
    }

    inline Bounds2 Expanded(Bounds2 bounds, float padding)
    {
        bounds.minX -= padding;
        bounds.minZ -= padding;
        bounds.maxX += padding;
        bounds.maxZ += padding;
        return bounds;
    }

    class Quadtree
    {
    public:
        struct Item
        {
            Bounds2 bounds;
            std::size_t index = 0;
        };

        explicit Quadtree(Bounds2 bounds, std::size_t capacity = 8, std::size_t maxDepth = 8)
            : root(std::make_unique<Node>(bounds, 0)),
              capacity(capacity),
              maxDepth(maxDepth)
        {
        }

        void Insert(Bounds2 bounds, std::size_t index)
        {
            Insert(*root, {bounds, index});
        }

        std::vector<std::size_t> Query(Bounds2 bounds) const
        {
            std::vector<std::size_t> result;
            Query(*root, bounds, result);
            return result;
        }

    private:
        struct Node
        {
            Node(Bounds2 bounds, std::size_t depth)
                : bounds(bounds), depth(depth)
            {
            }

            Bounds2 bounds;
            std::size_t depth = 0;
            std::vector<Item> items;
            std::array<std::unique_ptr<Node>, 4> children;
        };

        void Insert(Node& node, Item item)
        {
            if (node.children[0])
            {
                if (Node* child = ContainingChild(node, item.bounds))
                {
                    Insert(*child, item);
                    return;
                }
            }

            node.items.push_back(item);
            if (node.items.size() > capacity && node.depth < maxDepth)
                Subdivide(node);
        }

        void Query(const Node& node, Bounds2 bounds, std::vector<std::size_t>& result) const
        {
            if (!Intersects(node.bounds, bounds))
                return;

            for (const Item& item : node.items)
            {
                if (Intersects(item.bounds, bounds))
                    result.push_back(item.index);
            }

            for (const auto& child : node.children)
            {
                if (child)
                    Query(*child, bounds, result);
            }
        }

        void Subdivide(Node& node)
        {
            if (node.children[0])
                return;

            const float midX = (node.bounds.minX + node.bounds.maxX) * 0.5f;
            const float midZ = (node.bounds.minZ + node.bounds.maxZ) * 0.5f;
            node.children[0] = std::make_unique<Node>(Bounds2{node.bounds.minX, node.bounds.minZ, midX, midZ}, node.depth + 1);
            node.children[1] = std::make_unique<Node>(Bounds2{midX, node.bounds.minZ, node.bounds.maxX, midZ}, node.depth + 1);
            node.children[2] = std::make_unique<Node>(Bounds2{node.bounds.minX, midZ, midX, node.bounds.maxZ}, node.depth + 1);
            node.children[3] = std::make_unique<Node>(Bounds2{midX, midZ, node.bounds.maxX, node.bounds.maxZ}, node.depth + 1);

            std::vector<Item> retained;
            retained.reserve(node.items.size());
            for (const Item& item : node.items)
            {
                if (Node* child = ContainingChild(node, item.bounds))
                    Insert(*child, item);
                else
                    retained.push_back(item);
            }
            node.items = std::move(retained);
        }

        Node* ContainingChild(Node& node, Bounds2 bounds)
        {
            for (auto& child : node.children)
            {
                if (child && Contains(child->bounds, bounds))
                    return child.get();
            }
            return nullptr;
        }

        std::unique_ptr<Node> root;
        std::size_t capacity = 8;
        std::size_t maxDepth = 8;
    };
}
