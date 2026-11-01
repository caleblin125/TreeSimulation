#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <functional>
#include <vector>
#include <memory>

#include "simulation/physics.hpp"

// A limb is an elliptic cylinder along its local y axis. Width is the full
// cross-section extent along local x and thickness along local z, so the
// broad face of a flat limb faces local z.
class Limb : public std::enable_shared_from_this<Limb>
{
    typedef std::vector<std::shared_ptr<Limb>> LimbList;

public:
    typedef struct
    {
        glm::vec3 root;
        float length;
        float width;
        float thickness;
        glm::quat orientation;
        double ownWeight;
        double carriedWeight;
        glm::vec3 carriedTorque;
        double groundExtraction;
    } LimbData;

    Limb(glm::vec3 root, glm::quat orientation, float initialLength = 0.0f, float initialWidth = 0.0f, float initialThickness = 0.0f, double groundExtraction = 1.0);
    Limb(std::weak_ptr<Limb> parent, float attachment, glm::quat rotation, float initialLength = 0.0f, float initialWidth = 0.0f, float initialThickness = 0.0f, double groundExtraction = 1.0);
    void Spawn(float attachment, glm::quat localRotation, float initialLength, float initialWidth, float initialThickness);

    LimbData GetData()
    {
        return {root, length, width, thickness, orientation, ownWeight, carriedWeight, carriedTorque, groundExtraction};
    }

    std::vector<std::weak_ptr<Limb>> GetChildren()
    {
        std::vector<std::weak_ptr<Limb>> weak_list;
        weak_list.reserve(children.size());
        for (const auto &shared_ptr : children)
        {
            weak_list.push_back(shared_ptr);
        }
        return weak_list;
    }

    std::vector<std::shared_ptr<Limb>> CollectSubtree();
    void UpdatePose();
    // Gathers own weight, wind drag and everything the children carry. Returns carried weight.
    double UpdateLoads(double density, const Physics::Wind& wind = {});
    // Removes overloaded children and returns the weight that broke off.
    double PruneOverstressed(double bendingStrength);
    // Removes children for which `shouldShed` holds and returns the weight shed.
    double ShedWhere(const std::function<bool(const Limb&)>& shouldShed);
    bool ExceedsStrength(double bendingStrength) const;
    double BendingUtilization(double bendingStrength) const;
    double Volume() const;
    double SurfaceArea() const;
    double AboveGroundSurfaceArea() const;
    double UndergroundSurfaceArea() const;
    double GroundExtractionArea() const;
    // Shadow area of the above-ground part toward `direction`.
    double ProjectedArea(const glm::vec3& direction) const;
    double NutrientPathLength() const;
    unsigned int Depth() const;
    bool IsUnderground() const;
    bool CanPhotosynthesize() const;
    bool RayIntersects(const glm::vec3& origin, const glm::vec3& direction) const;
    void SetDimensions(float newLength, float newWidth, float newThickness);
    float GetLength() const;
    float GetWidth() const;
    float GetThickness() const;
    // Half the larger cross-section extent: how far the limb reaches off its axis.
    float GetReach() const;
    float GetAspectRatio() const;
    const glm::vec3& GetRoot() const;
    const glm::vec3& GetEnd() const;
    const glm::quat& GetOrientation() const;
    double GetCarriedWeight() const;
    const glm::vec3& GetCarriedTorque() const;
    const glm::vec3& GetOwnDrag() const;
    bool HasChildren() const;

    float ageYears = 0.0f;
    // Length when this limb last branched; it may branch again once it has grown past it.
    float lengthAtLastSplit = 0.0f;
    // Set once this segment has reached full length and opened a bud at its tip.
    bool hasLeader = false;
    // Inherited identity, 0 to 1, set by the parent when this limb budded and
    // fed to its growth brain, so a lineage of limbs can specialise.
    float branchClass = 0.0f;
    // Running average of the energy this limb earns minus what it costs, J per year.
    double energyBalance = 0.0;
    // Bending load relative to strength in a typical storm, worst direction.
    double designUtilization = 0.0;

protected:
    glm::vec3 root;   // Global root
    glm::vec3 end;    // Global end
    float attachment; // Position on parent [0, 1]

    glm::quat orientation; // Global orientation
    glm::quat rotation;    // Relative rotation

    float length = 0;
    float width = 0;
    float thickness = 0;

    std::weak_ptr<Limb> parent;
    LimbList children;
    double groundExtraction = 1.0;

    double ownWeight = 0.0;
    double carriedWeight = 0.0;
    glm::vec3 ownDrag = glm::vec3(0.0f);
    glm::vec3 carriedForce = glm::vec3(0.0f);
    glm::vec3 carriedTorque = glm::vec3(0.0f);

    friend class Organism;
};
